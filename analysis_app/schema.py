FIELDS = {'observation_id': 'str', 'machine': 'str', 'line': 'str', 'tick': 'int', 'residual': 'float', 'maintenance': 'int', 'visible_at': 'int', 'family': 'str'}
REQUIRED = ['observation_id', 'machine', 'line', 'residual']
ENTITY = 'machine'
GROUP = 'line'
MEASURE = 'residual'
PRODUCT = 'adaptive-fdr-monitor'
UNIT = '标准化残差'
