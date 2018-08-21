import os
from setuptools import setup, find_packages


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
        return fp.read()


setup(
    name='ACValidator',
    version="1.0.0",
    author="Shobana Sekar",
    author_email="ssekar@tgen.org",
    url="https://github.com/tgen/ACValidator",
    description="Assembly based approach for in silico validation of circular RNAs",
    keywords="",
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=2.7.13',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
    ],
    install_requires=[
        'pysam',
    ],
    #scripts=['ACV_launcher.sh'],
    entry_points={
        'console_scripts': [
            'ACValidator=ACValidator.ACValidator_v1:main',
        ],
    }
)
