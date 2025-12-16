import os

setup_config = os.environ['config']

if setup_config == 'CHEESEHEAD':
    from setup_CHEESEHEAD import*
elif setup_config == 'synthetic':
    from setup_synthetic import*
else:
    print('config undefined')
    
    
    
# Import other modules here