import setuptools

setuptools.setup(
    name="DrivAER",
    version="0.0.1",
    description="Identify driving transcriptional regulator in single cell embeddings",
    author='Lukas Simon, Fangfang Yan',
    author_email="Lukas.Simon@uth.tmc.edu",
    packages=['DrivAER'],
    install_requires=['sklearn','scanpy','anndata',
                      'pandas','seaborn','matplotlib','dca','tensorflow'
                      ],
    package_data={'DrivAER': ['data/*.txt','annotations/*.gmt','annotations/*.tsv']}
)