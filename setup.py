from setuptools import setup, find_packages

entry_points = """\
[console_scripts]
foo=bar:main
"""

setup(
    name='coolseq',
    version='0',
    author='Alexander Smith',
    author_email='asmitl@gmu.edu',
    url='https://github.com/ajsmith/sequence-analysis',
    packages=find_packages('src'),
    package_dir={'': 'src'},
#    entry_points=entry_points,
    install_requires=[
        'PyYAML',
    ],
)
