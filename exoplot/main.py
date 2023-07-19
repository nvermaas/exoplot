import os, sys, argparse
import json
from services.service_submit import do_submit
from services.service_processor import do_processor
from imaging.fits_imaging import draw_grid, draw_extra, get_min_max_ra_dec
from database.exoplanets import load_payload_from_database

def main():
    def get_arguments(parser):
        """
        Gets the arguments with which this application is called and returns
        the parsed arguments.
        If a argfile is give as argument, the arguments will be overrided
        The args.argfile need to be an absolute path!
        :param parser: the argument parser.
        :return: Returns the arguments.
        """
        args = parser.parse_args()
        if args.argfile:
            args_file = args.argfile
            if os.path.exists(args_file):
                parse_args_params = ['@' + args_file]
                # First add argument file
                # Now add command-line arguments to allow override of settings from file.
                for arg in sys.argv[1:]:  # Ignore first argument, since it is the path to the python script itself
                    parse_args_params.append(arg)
                print(parse_args_params)
                args = parser.parse_args(parse_args_params)
            else:
                raise (Exception("Can not find parameter file " + args_file))
        return args


    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')


    parser.add_argument("--astrometry_url",
                        default="http://nova.astrometry.net",
                        help="astrometry_api_key")
    parser.add_argument("--astrometry_api_key",
                        default=None,
                        help="astrometry_api_key")
    parser.add_argument("--source",
                        default="annotated",
                        help="draw on which image? 'raw', 'annotated'")
    parser.add_argument("--submission_id",
                        default=None,
                        help="if the submission_id is known then now new job is submitted, but instead the existing job is checked for results.")
    parser.add_argument("--title",
                        default="my title",
                        help="enter your own title for on the output image")
    parser.add_argument("--size",
                        default=20,
                        help="circle size")
    parser.add_argument("--path_to_image",
                        default="../examples/input_image.jpg",
                        help="filename of the input image")
    parser.add_argument("--output_dir",
                        default="../outputs",
                        help="directory where the output results will be stored")

    # All parameters in a file
    parser.add_argument('--argfile',
                        nargs='?',
                        type=str,
                        help='Ascii file with arguments (overrides all other arguments')

    args = get_arguments(parser)

    print(f"--- ExoPlot - 19 jul 2023 ---")


    # check for an apikey
    if not args.astrometry_api_key:
        print("An astrometry api_key is required to use this application. See https://nova.astrometry.net/api_help")
        exit(-1)

    # submit the image to astrometry.net
    submission_id = do_submit(args)

    # download the dataproducts from astrometry.net
    job_id = do_processor(submission_id, args)

    # draw exoplanets on the image
    path_to_fits_file = os.path.join(args.output_dir,f"{job_id}.fits")
    path, file = os.path.split(args.path_to_image)
    input_filename = os.path.splitext(file)

    # create <job>_grid.<ext>
    path_to_output_image = os.path.join(args.output_dir, f"{job_id}_grid{input_filename[1]}")
    draw_grid(path_to_fits_file, args.path_to_image, path_to_output_image, args.title, "degrees")

    # on which source file is the drawing done?
    if args.source == "annotated":
        source_image = os.path.join(args.output_dir,f"{job_id}_annotated{input_filename[1]}")
    elif args.source == "grid":
        source_image = os.path.join(args.output_dir,f"{job_id}_grid{input_filename[1]}")
    else:
        source_image = args.path_to_image

    # create <job>_exoplanets.<ext>
    path_to_output_image = os.path.join(args.output_dir, f"{job_id}_exoplanets{input_filename[1]}")

    # retrieve the sky coordinates based on the fits
    ra_start, ra_end, dec_start, dec_end = get_min_max_ra_dec(path_to_fits_file)

    # read the exoplanets from the database
    payload = load_payload_from_database(ra_start, ra_end, dec_start, dec_end, args)

    # draw the payload on the image
    draw_extra(path_to_fits_file, source_image, path_to_output_image, json.dumps(payload))


if __name__ == '__main__':
    main()
