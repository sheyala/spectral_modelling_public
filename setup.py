import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='spectral_modelling',
    version='0.0.1',
    author='Laura Cataldi',
    description='Spectral modelling tools',
    url='https://github.com/sheyala/spectral_modelling_public',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research'
    ],
    python_requires='==3.9.*,==3.10.*',
    install_requires=[
        'setuptools>=64',
        'numpy',
        'scipy==1.12.0',
        'shapely',
        'matplotlib',
        'pandas',
        'obspy==1.4.0'
    ]
)