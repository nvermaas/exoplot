# exoplot

Plotting exoplanets on your astronomy image (jpg or png).

### Description 
`exoplot` is a commandline tool that does 2 things
* plate solving: it uploads your image to http://nova.astrometry.net for 'plate-solving' and downloads the results. 
* plotting: it plots the contents of a database containing all confirmed exoplanets on your image.

**Important!** you need to request an api_key for the plate solving functionality: https://nova.astrometry.net/api_help


### github
https://github.com/nvermaas/exoplot

### pypi
https://test.pypi.org/project/exoplot/


### Installation

pip install -i https://test.pypi.org/simple/ exoplot==0.1.0

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
