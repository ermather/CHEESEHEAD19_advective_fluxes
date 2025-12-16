import os

setup_config = os.environ['config']

if setup_config == 'CHEESEHEAD':
    import setup_CHEESEHEAD
elif setup_config == 'synthetic':
    import setup_synthetic
else:
    print('config undefined')