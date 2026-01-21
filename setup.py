from setuptools import setup, find_packages
import versioneer


def readme():
    with open('README.md', 'r') as ip:
        return ip.read()


setup(
    name='recon-pop',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=('Recon: Reconstruction of Estimated Communities from '
                 'Observed Numbers'),
    long_description=readme(),
    url='https://github.com/ArnaoutLab/Recon',
    author='Arnaout Lab',
    license='Apache License, Version 2.0',
    classifiers=['Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3'],
    packages=find_packages(),
    install_requires=['numpy', 'scipy'],
    entry_points={'console_scripts': ['Recon = Recon.cli:main']})
