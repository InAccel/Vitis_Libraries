"""Setup inaccel-vitis package."""
from setuptools import find_namespace_packages, setup

setup(
    name = 'inaccel-vitis',
    packages = find_namespace_packages(include = ['inaccel.*']),
    namespace_packages = ['inaccel'],
    version = '0.2',
    license = 'Apache-2.0',
    description = 'InAccel Vitis Libraries',
    author = 'InAccel',
    author_email = 'info@inaccel.com',
    url = 'https://docs.inaccel.com',
    keywords = ['InAccel Coral', 'FPGA', 'inaccel', 'Vitis'],
    install_requires = [
        'coral-api==2.*', 'opencv-python',
    ],
    include_package_data = True,
    zip_safe = True,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires = '>=3.7',
)
