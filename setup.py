from setuptools import setup, find_packages


setup(
    name='proxysiphon',
    version='0.0.1a0',
    description='Internal lab tool to parse and clean marine sediment proxy data.',
    license='GPLv3',

    author='S. Brewster Malevich',
    author_email='malevich@email.arizona.edu',
    url='https://github.com/brews/proxysiphon',
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
    ],
    keywords='marine paleoclimate',

    packages=find_packages(exclude=['docs']),

    install_requires=['numpy', 'scipy', 'matplotlib', 'pandas', 'cartopy', 'attrs',
                      'carbonferret', 'snakebacon', 'erebus'],
    tests_require=['pytest'],
    package_data={'proxysiphon': ['tests/*.txt']},
)
