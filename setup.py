from setuptools import setup


setup(
    # Self-descriptive entries which should always be present
    name='cgnp',
    author='Caroline',
    author_email='cjs323@lehigh.edu',
    license='MIT',
    version='0.0.0',
    description='Builds a silica, alkane-coated, coarse-grained nanoparticle.',
    zip_safe=False,
    entry_points={
        'mbuild.plugins':[
        "cgnp = cgnp.cgnp:cgnp"
        ]
        }
    )
