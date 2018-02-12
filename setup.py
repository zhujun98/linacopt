#!/usr/bin/env python
from setuptools import setup, find_packages


if __name__ == '__main__':
    setup(
        name='linacopt',
        version='1.0.0',
        packages=find_packages(),

        install_requires=[
            'numpy>=1.13',
            'pandas>=0.20.3',
            'matplotlib>=2.1.0',
            'scipy>=1.00'
        ],

        author='Jun Zhu',
        author_email='zhujun981661@gmail.com',
        url='https://github.com/zhujun98/linacopt',
        download_url='https://github.com/zhujun98/linacopt',
        description='Python API for beam dynamics optimization',
        long_description='linacopt a Python API for local and '
                         'global optimization of linac beam dynamics '
                         'simulation with ASTRA and IMPACT-T',
        license='GNU',

        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        classifiers=['Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: Developers',
                     'License :: GNU',
                     'Programming Language :: Python :: 3.6',
                     'Topic :: Scientific/Engineering',
                     'Topic :: Software Development',
                     'Operating System :: Microsoft :: Windows',
                     'Operating System :: POSIX',
                     'Operating System :: Unix',
                     'Operating System :: MacOS',
                     ]
    )
