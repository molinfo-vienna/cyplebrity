from setuptools import find_packages, setup

setup(
    name="cyplebrity",
    version="0.3.2",
    maintainer="Johannes Kirchmair",
    maintainer_email="johannes.kirchmair@univie.ac.at",
    packages=find_packages(),
    url="https://github.com/molinfo-vienna/cyplebrity",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="BSD 3-Clause License",
    include_package_data=True,
    install_requires=[
        "rdkit==2020.09.1",
        "scikit_learn==0.23.2",
        "numpy==1.19.2",
        "molvs==0.1.1",
        "nerdd-module>=0.3.16",
        "fpsim2==0.2.8",
        # avoid warnings about numpy.distutils
        "setuptools < 60.0",
        # install importlib-resources and importlib-metadata for old Python versions
        "importlib-resources>=5; python_version<'3.9'",
        "importlib-metadata>=4.6; python_version<'3.10'",
    ],
    extras_require={
        "dev": ["mypy", "ruff"],
        "test": [
            "pytest",
            "pytest-watch",
            "pytest-cov",
            "pytest-bdd==7.3.0",
            "hypothesis",
            "hypothesis-rdkit",
        ],
    },
    entry_points={
        "console_scripts": [
            "cyplebrity=cyplebrity.__main__:main",
        ],
    },
)
