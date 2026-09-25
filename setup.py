from setuptools import setup, find_packages
setup(name="digitalembryo", version="0.1", package_dir={"": "src"},
      packages=find_packages("src"),
    entry_points={'console_scripts': ['embryosim=digitalembryo.cli:main']},
)
