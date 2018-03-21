from setuptools import setup, find_packages

setup(
    name = "cpssdk",
    version = "0.0.1",
    keywords = ("pip", "datacanvas", "eds", "xiaoh"),
    description = "cps time processing package",
    long_description = "read time and calculate trip duration package in cyber physical system processing",
    license = "MIT Licence",
    
    url = "http://zhihanfang.com",
    author = "ogre",
    author_email = "ezhihan@gmail.com",
    
    packages = find_packages(),
    include_package_data = True,
    platforms = "any",
    install_requires = ['datetime','simplejson']
)
