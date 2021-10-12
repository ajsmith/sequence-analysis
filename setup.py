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
    packages=['coolseq'],
    package_dir={'': 'src'},
#    entry_points=entry_points,
    include_package_data=True,
    package_data={
        'coolseq': ['align/samples.fasta'],
    },
    install_requires=[
        'PyYAML',
        'numpy',
    ],
)
