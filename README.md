# exoplot

Plotting exoplanets on your astronomy image (jpg or png).
This is a commandline tool in Python that does 2 things:
* plate solving: it uploads your image to http://nova.astrometry.net and downloads the results. 
* plotting: it plots the contents of an exoplanet database on a copy of your image.

**Important!** you need to request an api_key for the plate solving functionality: https://nova.astrometry.net/api_help

### Resources
* github (source code): https://github.com/nvermaas/exoplot
* pypi (pip installable): https://test.pypi.org/project/exoplot/

### Install
```
> virtualenv env
> source env/bin/activate
> pip install -i https://test.pypi.org/simple/ exoplot
```

### Usage
> exoplot -h

----

## For Developers

### clone the repo
> clone https://github.com/nvermaas/exoplot.git

### run locally (like in PyCharm)
* Execute the ```main.py``` file
* Required parameter: ```--astrometry_api_key <your-api-key>
* See the options: `-h`

### build
> py -m build
### upload to pypi
> py -m twine upload --repository testpypi dist/*

### deploy
> pip install -i https://test.pypi.org/simple/ exoplot
