from setuptools import setup, find_packages


setup(
    name='proxysiphon',
    version='0.0.1a1',
    description='Internal lab tool to parse and clean marine sediment proxy data.',
    license='GPLv3',

    author='S. Brewster Malevich',
    author_email='malevich@email.arizona.edu',
    url='https://github.com/brews/proxysiphon',
    classifiers=[
        'Development Status :: 1 - Planning',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        'Programming Language :: Python :: 3',
    ],
    keywords='marine paleoclimate',

    packages=find_packages(exclude=['docs']),

    install_requires=['numpy', 'scipy', 'matplotlib', 'pandas', 'cartopy',
                      'chardet', 'carbonferret', 'snakebacon', 'erebusfall'],
    tests_require=['pytest'],
    package_data={'proxysiphon': ['tests/*.txt']},
)
