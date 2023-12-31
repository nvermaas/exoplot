# exoplot

Plotting exoplanets on your astronomy image
* plate solving by uploading your image to http://nova.astrometry.net and downloads the results. 
* plotting all confirmed exoplanets (as of july 2023) on a copy of your image.

**Important!** 
you need to request an api_key for the plate solving functionality: https://nova.astrometry.net/api_help

### Resources
* github (source code): https://github.com/nvermaas/exoplot
* pypi (pip installable): https://pypi.org/project/exoplot


### Install
```
> virtualenv env
> source env/bin/activate
> pip install astropy Pillow
> pip install exoplot --upgrade
```

### Usage
> exoplot -h

```
usage: main.py [-h] [--astrometry_api_key ASTROMETRY_API_KEY]
               [--astrometry_url ASTROMETRY_URL]
               [--exoplanets_db EXOPLANETS_DB] [--source SOURCE]
               [--submission_id SUBMISSION_ID] [--title TITLE] [--size SIZE]
               [--path_to_image PATH_TO_IMAGE] [--output_dir OUTPUT_DIR]
               [--argfile [ARGFILE]]

options:
  -h, --help            show this help message and exit
  --astrometry_api_key ASTROMETRY_API_KEY
                        the api_key that you requested at
                        https://nova.astrometry.net/api_help
  --astrometry_url ASTROMETRY_URL
                        astrometry instance for plate solving. If not provided
                        http://nova.astrometry.net will be used
  --exoplanets_db EXOPLANETS_DB
                        path to exoplanets sqlite3 database. If not provided
                        the onboard (july 2023) database will be used
  --source SOURCE       Background image to draw exoplanets on? Options are:
                        'raw', 'annotated', 'grid'
  --submission_id SUBMISSION_ID
                        if the submission_id is known then no new job is
                        submitted, but instead the existing job is checked for
                        results. Useful for playing with the other parameters
                        while keeping the existing plate solving results.
  --title TITLE         display your title on the output image
  --size SIZE           circle size in pixels
  --path_to_image PATH_TO_IMAGE
                        input image, if not provided the onboard example image
                        will be used
  --output_dir OUTPUT_DIR
                        directory where the output results will be stored
  --argfile [ARGFILE]   Ascii file with arguments (overrides all other
                        arguments

```

### The most important options

#### --astrometry_api_key
you need to request an api_key for the plate solving functionality: https://nova.astrometry.net/api_help

#### --path_to_image
This is the path to your .jpg or .png file.
If not provided, it will use an example image of Lyra

#### --output_dir
The directory where the results of the plate solving and image drawing are stored.
All files start with a job_id. The following files will be there after a run.
(note: for size considerations, the following images are partial screenshots)

* <job_id>.fits : the original input image converted to fits
* <job_id>_annotated.jpg : constellation lines, stars and dso names, downloaded from astrometry.net
![](/docs/annotated.jpg)
* <job_id>_grid.jpg : the previous annotated image, but with a coordinate grid drawn (by the exoplot tool) 
![](/docs/grid.jpg)
* <job_id>_sky_globe.jpg : a global skymap showing the location of your image 
![](/docs/sky_globe.jpg)
* <job_id>_sky_plot.jpg : a more detailed skymap showing the location of your image
![](/docs/sky_plot.jpg)
* <job_id>_exoplanets.jpg : a copy of your image with all currently known confirmed exoplanets plotted.
![](/docs/exoplanets.jpg)

#### --source
You can choose on which image you plot the exoplanets:
* your **raw** image
* the plate solved **annotated** image
* the annotated image with a RA/dec **grid** drawn on it.


----

## For Developers

### clone the repo
> clone https://github.com/nvermaas/exoplot.git

### run locally (like in PyCharm)
* Execute the ```main.py``` file
* Required parameter: ```--astrometry_api_key <your-api-key>```
* See the options: `-h`

### build
> py -m build
### upload to pypi
> py -m twine upload dist/*

### deploy
> pip install exoplot --upgrade


*disclaimer*:
I originally used an old version of an astrometry.net client, but I cannot remember where I found it.
So the following class was not written by me, but by some anonymous hero.
https://github.com/nvermaas/exoplot/blob/master/exoplot/astrometry/astrometry_client.py
