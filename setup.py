from setuptools import setup, find_packages

with open('README.rst') as f:
    long_description = ''.join(f.readlines())


setup(
    name='pairef',
    version='1.2.1',
    description='Automatic PAIRed REFinement protocol',
    long_description=long_description,
    author='Martin Maly',
    author_email='martin.maly@fjfi.cvut.cz',
    keywords='macromolecular crystallography, research',
    license='GNU Lesser General Public License v3 (LGPLv3)',
    url='https://pairef.fjfi.cvut.cz/',
    packages=find_packages(),
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'pairef = pairef.launcher:run_pairef',
            'pairef-gui = pairef.gui:gui',
        ]
    },
    install_requires=['numpy', 'matplotlib'],
    package_data={'pairef': ['static/*.css']},
    data_files=[('bitmaps', ['pairef/static/pairef_logo_64.png'])],
    include_package_data=True,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
