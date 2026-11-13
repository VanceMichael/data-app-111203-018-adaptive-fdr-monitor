import copy
import csv
import json
import math
import random
import tempfile
import unittest
from pathlib import Path

from analysis_app.fdr import (GAMMA_C, LORD_INITIAL_FRACTION,
                              OnlineFDRController, Stream, gamma_at)
from analysis_app.io import FIELDS, load_records
from analysis_app.monitor import (Monitor, MonitorConfig,
                                  DISPOSITION_LATE, DISPOSITION_DUPLICATE)
from analysis_app.references import ReferenceBook
from analysis_app.stats import lag1_autocorrelation, t_two_sided_p


def row(observation_id, machine='M00', family='F0', tick=0, residual=0.0,
        maintenance=0, visible_at=None, line='L0'):
    return {'observation_id': observation_id, 'machine': machine,
            'family': family, 'line': line, 'tick': tick,
            'residual': residual, 'maintenance': maintenance,
            'visible_at': tick if visible_at is None else visible_at}


class StatsTest(unittest.TestCase):
    def test_t_table_critical_values(self):
        self.assertAlmostEqual(t_two_sided_p(2.228, 10), 0.05, places=4)
        self.assertAlmostEqual(t_two_sided_p(3.169, 10), 0.01, places=4)
        self.assertAlmostEqual(t_two_sided_p(2.776, 4), 0.05, places=4)
        self.assertEqual(t_two_sided_p(0.0, 5), 1.0)

    def test_t_p_monotone_and_bounded(self):
        previous = 1.0
        for t in [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]:
            p = t_two_sided_p(t, 8)
            self.assertLess(p, previous)
            self.assertGreater(p, 0.0)
            previous = p

    def test_lag1_autocorrelation(self):
        self.assertIsNone(lag1_autocorrelation([1.0]))
        self.assertEqual(lag1_autocorrelation([2.0, 2.0, 2.0]), 0.0)
        self.assertAlmostEqual(lag1_autocorrelation([1.0, 2.0, 3.0, 4.0]), 0.25)


class FDRControllerTest(unittest.TestCase):
    def test_gamma_summable(self):
        total = math.fsum(gamma_at(t) for t in range(1, 100000))
        self.assertAlmostEqual(total, 1.0, places=4)
        self.assertAlmostEqual(GAMMA_C, 6.0 / math.pi ** 2)

    def test_lord3_first_threshold_is_initial_wealth(self):
        controller = OnlineFDRController(0.05, Stream.LORD3)
        self.assertAlmostEqual(controller.threshold(),
                               LORD_INITIAL_FRACTION * 0.05)

    def test_lord3_wealth_replenishes_on_rejection(self):
        controller = OnlineFDRController(0.05, Stream.LORD3)
        first = controller.submit(0.0001)
        self.assertTrue(first['rejected'])
        self.assertGreater(controller.wealth,
                           LORD_INITIAL_FRACTION * 0.05)
        second = controller.submit(1.0)
        self.assertFalse(second['rejected'])
        self.assertLess(controller.wealth,
                        first['wealth_after'])

    def test_lord3_never_recovers_from_exact_zero_without_floor(self):
        controller = OnlineFDRController(0.05, Stream.LORD3)
        for _ in range(200):
            controller.submit(1.0)
        self.assertGreater(controller.threshold(), 0.0)

    def test_bonferroni_thresholds_summable_bound(self):
        controller = OnlineFDRController(0.05, Stream.ONLINE_BONFERRONI)
        spent = math.fsum(controller.submit(0.0)['alpha'] for _ in range(5000))
        # 并集界：所有阈值之和不超过 envelope alpha
        self.assertLessEqual(spent, 0.05 + 1e-12)

    def test_bonferroni_rejects_under_arbitrary_dependence_schedule(self):
        controller = OnlineFDRController(0.1, Stream.ONLINE_BONFERRONI)
        entry = controller.submit(1e-9)
        self.assertTrue(entry['rejected'])
        self.assertIsNone(entry['wealth_after'])

    def test_invalid_p_value_rejected(self):
        controller = OnlineFDRController(0.05, Stream.LORD3)
        with self.assertRaises(ValueError):
            controller.submit(1.2)

    def test_restore_matches_running_state(self):
        controller = OnlineFDRController(0.05, Stream.LORD3)
        for p in [0.001, 0.9, 0.2, 0.0001, 0.8]:
            controller.submit(p)
        restored = OnlineFDRController.restore(controller.status().to_dict())
        self.assertEqual(restored.position, controller.position)
        self.assertEqual(restored.wealth, controller.wealth)
        self.assertEqual(restored.reject_positions, controller.reject_positions)
        fresh = OnlineFDRController(0.05, Stream.LORD3)
        for p in [0.001, 0.9, 0.2, 0.0001, 0.8, 0.3]:
            fresh.submit(p)
        restored.submit(0.3)
        self.assertEqual(restored.wealth, fresh.wealth)


class ReferenceBookTest(unittest.TestCase):
    def test_windows_never_mix_across_maintenance(self):
        book = ReferenceBook(reference_size=3)
        for tick, residual in enumerate([0.0, 0.1, -0.1, 0.2, 5.0]):
            book.observe(row(f'a{tick}', tick=tick, residual=residual))
        book.observe(row('maint', tick=5, maintenance=1))
        for tick, residual in [(6, 2.0), (7, 2.1), (8, 1.9)]:
            epoch, status = book.observe(row(f'b{tick}', tick=tick, residual=residual))
        ref0 = book.reference_for('M00', 0)
        ref1 = book.reference_for('M00', 1)
        self.assertEqual(ref0.start_tick, 0)
        self.assertEqual(ref1.start_tick, 6)
        self.assertAlmostEqual(ref1.mean, 2.0, places=6)
        self.assertNotIn(5.0, ref1.observation_ids)

    def test_cold_start_when_window_closes_early(self):
        book = ReferenceBook(reference_size=4)
        for tick in range(3):
            book.observe(row(f'w{tick}', tick=tick, residual=float(tick)))
        epoch, status = book.observe(row('m', tick=3, maintenance=1))
        self.assertIsNone(book.reference_for('M00', epoch))
        self.assertTrue(book.epoch_state('M00', epoch).closed)

    def test_zero_variance_window_does_not_freeze(self):
        book = ReferenceBook(reference_size=3)
        for tick in range(3):
            book.observe(row(f'z{tick}', tick=tick, residual=1.5))
        self.assertIsNone(book.reference_for('M00', 0))

    def test_snapshot_roundtrip(self):
        book = ReferenceBook(reference_size=3)
        for tick, residual in enumerate([0.0, 0.1, -0.1, 0.2]):
            book.observe(row(f's{tick}', tick=tick, residual=residual))
        restored = ReferenceBook.restore(book.snapshot())
        self.assertEqual(restored.reference_for('M00', 0).mean,
                         book.reference_for('M00', 0).mean)


class MonitorOrderingTest(unittest.TestCase):
    def setUp(self):
        self.config = MonitorConfig(reference_size=3,
                                    family_weights={'F0': 0.5, 'F1': 0.5})

    def test_sorted_by_visible_at_regardless_of_input_order(self):
        rows = [row('late', tick=5, visible_at=5, residual=2.0),
                row('early', tick=0, visible_at=0, residual=0.1)]
        monitor = Monitor(self.config)
        monitor.ingest(rows)
        self.assertEqual([d['observation_id'] for d in monitor.dispositions],
                         ['early', 'late'])

    def test_late_arrival_quarantined_never_spends_budget(self):
        monitor = Monitor(self.config)
        monitor.ingest([row('a', tick=10, visible_at=10, residual=1.0),
                        row('b', tick=11, visible_at=11, residual=1.0)])
        before = monitor.global_position
        monitor.ingest([row('c', tick=1, visible_at=2)])
        self.assertEqual(monitor.global_position, before)
        self.assertEqual(monitor.dispositions[-1]['status'],
                         DISPOSITION_LATE)

    def test_pending_buffer_orders_within_batch_before_watermark(self):
        # 同批先给晚可见、再给早可见：缓冲排序后早的先处理，不产生迟到
        monitor = Monitor(MonitorConfig(
            reference_size=2, family_weights={'F0': 1.0}))
        batch = [row('late-ish', tick=3, visible_at=3, residual=0.9),
                 row('early-1', tick=1, visible_at=1, residual=0.2),
                 row('early-0', tick=0, visible_at=0, residual=0.1),
                 row('mid', tick=2, visible_at=2, residual=-0.3)]
        monitor.ingest(batch)
        statuses = [d['status'] for d in monitor.dispositions]
        self.assertNotIn(DISPOSITION_LATE, statuses)
        self.assertEqual(
            [d['visible_at'] for d in monitor.dispositions], [0, 1, 2, 3])

    def test_watermark_prescreen_isolates_late_in_mixed_batch(self):
        monitor = Monitor(MonitorConfig(
            reference_size=2, family_weights={'F0': 1.0}))
        monitor.ingest([row('seen', tick=5, visible_at=5, residual=0.5)])
        monitor.ingest([row('new-high', tick=6, visible_at=6, residual=0.6),
                        row('new-low', tick=0, visible_at=0, residual=0.0)])
        statuses = {d['observation_id']: d['status']
                    for d in monitor.dispositions}
        self.assertEqual(statuses['new-low'], DISPOSITION_LATE)
        self.assertNotEqual(statuses['new-high'], DISPOSITION_LATE)

    def test_duplicate_does_not_reconsume_budget(self):
        monitor = Monitor(self.config)
        monitor.ingest([row('a', tick=0, residual=0.0),
                        row('b', tick=1, residual=0.1),
                        row('c', tick=2, residual=-0.1),
                        row('d', tick=3, residual=0.5)])
        ledger_len = len(monitor.budget_ledger)
        monitor.ingest([row('d', tick=3, residual=0.5)])
        self.assertEqual(len(monitor.budget_ledger), ledger_len)
        self.assertEqual(monitor.dispositions[-1]['status'],
                         DISPOSITION_DUPLICATE)

    def test_batched_ingest_matches_one_shot(self):
        records = load_records('data/observations.csv')
        ordered = sorted(records, key=Monitor._order_key)
        one_shot = Monitor(self.config)
        one_shot.ingest(records)
        batched = Monitor(self.config)
        half = len(ordered) // 2
        batched.ingest(ordered[:half])
        batched = Monitor.restore(batched.state_snapshot())
        batched.ingest(ordered[half:])
        self.assertEqual(batched.budget_ledger, one_shot.budget_ledger)
        self.assertEqual(batched.watermark, one_shot.watermark)

    def test_global_positions_are_contiguous_and_immutable(self):
        records = load_records('data/observations.csv')
        monitor = Monitor(self.config)
        monitor.ingest(records)
        positions = [entry['global_position'] for entry in monitor.budget_ledger]
        self.assertEqual(positions, list(range(1, len(positions) + 1)))
        snapshot = copy.deepcopy(monitor.state_snapshot())
        monitor.ingest(records)  # all duplicates
        self.assertEqual(
            [entry['global_position'] for entry in monitor.budget_ledger],
            [entry['global_position'] for entry in snapshot['budget_ledger']])


class MonitorReferenceTest(unittest.TestCase):
    def test_maintenance_keeps_global_ledger_and_rebuilds_reference(self):
        records = load_records('data/observations.csv')
        config = MonitorConfig(reference_size=3,
                               family_weights={'F0': 0.5, 'F1': 0.5})
        monitor = Monitor(config)
        monitor.ingest(records)
        ledger_before = len(monitor.budget_ledger)
        refs = monitor.reference_book.snapshot()['machines']['M00']
        self.assertEqual([item['epoch'] for item in refs], [0, 1])
        self.assertLessEqual(refs[0]['reference']['end_tick'], 5)
        self.assertEqual(refs[1]['reference']['start_tick'], 7)
        post = [entry for entry in monitor.budget_ledger
                if entry['machine'] == 'M00' and entry['epoch'] == 1]
        self.assertTrue(post)
        self.assertTrue(all(
            entry['validity_conditions']['reference_start_tick'] == 7
            for entry in post))
        self.assertEqual(len(monitor.budget_ledger), ledger_before)

    def test_cold_start_state_when_reference_insufficient(self):
        records = load_records('data/observations.csv')
        config = MonitorConfig(reference_size=10,
                               family_weights={'F0': 0.5, 'F1': 0.5})
        monitor = Monitor(config)
        monitor.ingest(records)
        monitor.finalize()
        self.assertEqual(monitor.global_position, 0)
        cold = [d for d in monitor.dispositions if d['status'] == 'cold_start']
        self.assertEqual(len(cold), 264)

    def test_unknown_family_is_quarantined(self):
        records = load_records('data/observations.csv')
        monitor = Monitor(MonitorConfig(family_weights={'F0': 1.0}))
        monitor.ingest(records)
        self.assertTrue(all(
            d['status'] == 'quarantined' for d in monitor.dispositions
            if d['family'] == 'F1'))

    def test_force_conservative_uses_only_bonferroni(self):
        records = load_records('data/observations.csv')
        config = MonitorConfig(reference_size=3, force_conservative=True,
                               family_weights={'F0': 0.5, 'F1': 0.5})
        monitor = Monitor(config)
        monitor.ingest(records)
        self.assertTrue(monitor.budget_ledger)
        self.assertTrue(all(entry['stream'] == 'online_bonferroni'
                            for entry in monitor.budget_ledger))


class MonitorRecoveryTest(unittest.TestCase):
    def test_state_roundtrip_preserves_everything(self):
        records = load_records('data/observations.csv')
        config = MonitorConfig(reference_size=3,
                               family_weights={'F0': 0.5, 'F1': 0.5})
        monitor = Monitor(config)
        monitor.ingest(records[:len(records) // 2])
        monitor.finalize()
        restored = Monitor.restore(monitor.state_snapshot())
        monitor.ingest(records[len(records) // 2:])
        monitor.finalize()
        restored.ingest(records[len(records) // 2:])
        restored.finalize()
        self.assertEqual(restored.budget_ledger, monitor.budget_ledger)
        self.assertEqual(restored.discoveries, monitor.discoveries)
        self.assertEqual(restored.watermark, monitor.watermark)

    def test_invalid_restored_config_rejected(self):
        monitor = Monitor(MonitorConfig(global_alpha=0.05))
        payload = monitor.state_snapshot()
        payload['config']['global_alpha'] = 1.5
        with self.assertRaises(ValueError):
            Monitor.restore(payload)
        payload = monitor.state_snapshot()
        payload['config']['lord_share'] = 2.0
        with self.assertRaises(ValueError):
            Monitor.restore(payload)


class HierarchicalBudgetTest(unittest.TestCase):
    def test_family_envelopes_follow_fixed_weights(self):
        config = MonitorConfig(global_alpha=0.08,
                               family_weights={'F0': 0.75, 'F1': 0.25})
        monitor = Monitor(config)
        monitor._ensure_family('F0')
        monitor._ensure_family('F1')
        lord0 = monitor.families['F0'].controllers['lord3']
        lord1 = monitor.families['F1'].controllers['lord3']
        self.assertAlmostEqual(lord0.envelope_alpha, 0.08 * 0.75 * 0.8)
        self.assertAlmostEqual(lord1.envelope_alpha, 0.08 * 0.25 * 0.8)
        # 两个 family 流包络之和不超过全局 alpha
        total = math.fsum(
            controller.envelope_alpha
            for runtime in monitor.families.values()
            for controller in runtime.controllers.values())
        self.assertLessEqual(total, 0.08 + 1e-12)

    def test_ledger_records_p_value_and_validity_conditions(self):
        records = load_records('data/observations.csv')
        monitor = Monitor(MonitorConfig(
            reference_size=3, family_weights={'F0': 0.5, 'F1': 0.5}))
        monitor.ingest(records)
        entry = monitor.budget_ledger[0]
        self.assertIn('p_value', entry)
        self.assertIn('null_hypothesis', entry)
        conditions = entry['validity_conditions']
        self.assertEqual(conditions['reference_size'], 3)
        self.assertEqual(len(conditions['reference_observation_ids']), 3)
        self.assertIn('stream_assumptions', conditions)
        self.assertFalse(conditions['degraded'])

    def test_equal_weights_freeze_and_late_family_quarantine(self):
        records = [row(f'a{i}', family='F0', tick=i, residual=0.1 * i)
                   for i in range(4)]
        monitor = Monitor(MonitorConfig(reference_size=3))
        monitor.ingest(records)
        self.assertEqual(monitor.families['F0'].weight, 1.0)
        restored = Monitor.restore(monitor.state_snapshot())
        restored.ingest([row('b1', machine='M1', family='F9', tick=0,
                             visible_at=4)])
        self.assertEqual(restored.dispositions[-1]['status'], 'quarantined')

    def test_equal_split_within_first_batch_sums_to_global_alpha(self):
        records_f0 = [row(f'a{i}', family='F0', tick=i, residual=0.1 * i)
                      for i in range(4)]
        records_f1 = [dict(item, observation_id=item['observation_id'] + 'b',
                           family='F1') for item in records_f0]
        monitor = Monitor(MonitorConfig(reference_size=3))
        monitor.ingest(records_f0 + records_f1)
        self.assertAlmostEqual(monitor.families['F0'].weight, 0.5)
        self.assertAlmostEqual(monitor.families['F1'].weight, 0.5)
        total = math.fsum(
            controller.envelope_alpha
            for runtime in monitor.families.values()
            for controller in runtime.controllers.values())
        self.assertAlmostEqual(total, 0.05)


class LongRunFDRTest(unittest.TestCase):
    """独立零假设下的长期误报表现（固定种子，确定性）。"""

    def test_iid_null_false_rate_bounded(self):
        rng = random.Random(99)
        false_alarms = 0
        tests = 0
        for trial in range(6):
            rows = []
            for tick in range(300):
                for m in range(8):
                    rows.append(row(f't{trial}-m{m}-{tick}',
                                    machine=f'M{m}', family=f'F{m % 2}',
                                    tick=tick, residual=rng.gauss(0, 1),
                                    visible_at=tick))
            monitor = Monitor(MonitorConfig(
                reference_size=20, ar_min_samples=60,
                family_weights={'F0': 0.5, 'F1': 0.5}))
            monitor.ingest(rows)
            monitor.finalize()
            false_alarms += len(monitor.discoveries)
            tests += monitor.global_position
        self.assertGreater(tests, 10000)
        # 独立零假设下 FDR 控制在 0.05；实际保守得多
        self.assertLess(false_alarms / tests, 0.05)

    def test_ar1_process_routes_to_bonferroni(self):
        rng = random.Random(7)
        phi = 0.6
        state = {m: 0.0 for m in range(6)}
        rows = []
        for tick in range(300):
            for m in range(6):
                noise = rng.gauss(0, math.sqrt(1 - phi ** 2))
                state[m] = phi * state[m] + noise
                rows.append(row(f'ar-m{m}-{tick}', machine=f'M{m}',
                                family=f'F{m % 2}', tick=tick,
                                residual=state[m], visible_at=tick))
        monitor = Monitor(MonitorConfig(
            reference_size=20, ar_min_samples=60, ar_cutoff=0.5,
            family_weights={'F0': 0.5, 'F1': 0.5}))
        monitor.ingest(rows)
        monitor.finalize()
        degraded = monitor.current_state()['degraded_machine_epochs']
        self.assertTrue(degraded)
        bonferroni = sum(1 for entry in monitor.budget_ledger
                         if entry['stream'] == 'online_bonferroni')
        self.assertGreater(bonferroni, 0)
        # 相关零假设下误报仍受全局 alpha 约束
        self.assertLessEqual(len(monitor.discoveries),
                             max(1, int(0.05 * monitor.global_position) + 2))


if __name__ == '__main__':
    unittest.main()
