import math
import statistics
from collections import Counter


def describe(records):
    return {'absolute_residual_sum': math.fsum(abs(row['residual']) for row in records), 'maintenance_records': sum(row['maintenance'] for row in records), 'machines': len({row['machine'] for row in records})}
