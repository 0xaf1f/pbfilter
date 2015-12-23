from setuptools import setup, find_packages

setup(
    name = 'pbfilter',
    version = '1.2',
    author = 'pbiDevNet',
    author_email = 'pbiDevNet@pacificbiosciences.com',
    license = 'LICENSE.txt',
    packages = find_packages('.'),
    package_dir = {'':'.'},
    zip_safe = False,
    entry_points = dict(console_scripts=[
        'filter_artifacts = pbfilter.filter_artifacts:main',
        'filter_plsh5 = pbfilter.filter_plsh5:main'
    ]),
    install_requires = [
        'numpy',
        'h5py',
        'pbcore'
    ]
)
