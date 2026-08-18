from setuptools import setup, find_packages

setup(
    name="zbox_nn_engine",
    version="1.0.0",
    author="Thijs Van der Snickt et al.",
    description="Topology-aware thermodynamics for DNA probe design",
    long_description=(
        "A standalone Python package implementing three scoring layers "
        "for DNA probe selectivity prediction: S_ECI (uniform-weight "
        "box-counting), N_box (analytical sub-box count), and Z_box,NN "
        "(NN-weighted partition function). Based on SantaLucia 2004 "
        "unified nearest-neighbor parameters."
    ),
    long_description_content_type="text/plain",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
