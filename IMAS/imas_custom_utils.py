import argparse
import imas

def parse_and_open_imas_entry():
    parser = argparse.ArgumentParser(description='---- Display scenario')
    parser.add_argument('-s', '--shot', help='Shot number', required=True, type=int)
    parser.add_argument('-r', '--run', help='Run number', required=True, type=int)
    parser.add_argument('-u', '--user_or_path', help='User or absolute path name where the data-entry is located', required=False, default='public', type=str)
    parser.add_argument('-d', '--database', help='Database name where the data-entry is located', required=False, default='iter', type=str)
    parser.add_argument('-b', '--backend', help='Backend (HDF5, MDSPLUS)', required=False, default='HDF5', type=str)
    parser.add_argument('-v', '--version', help='Version of data dictionary (3,4,...)', required=False, default='4', type=str)
    parser.add_argument('-q', '--quantity', help='Quantity to plot (Ip/BR/Jtor/Te)', required=False, default='4', type=str)
    parser.add_argument('-it', '--time_index', help='Time slice index', required=False, default=0, type=int)

    args = vars(parser.parse_args())

    shot = args["shot"]
    run = args["run"]
    user = args['user_or_path']
    database = args['database']
    dd_version = args['version']

    # Determine backend
    if args['backend'] == 'MDSPLUS':
        backend_id = imas.imasdef.MDSPLUS_BACKEND
    else:
        backend_id = imas.imasdef.HDF5_BACKEND

    uri = imas.DBEntry.build_uri_from_legacy_parameters(
        backend_id=backend_id,
        pulse=shot,
        run=run,
        db_name=database,
        user_name=user,
        data_version=dd_version
    )

    imas_entry = imas.DBEntry(uri=uri, mode='r')
    imas_entry.open()
    
    output = {}
    output["entry"] = imas_entry
    output["args"] = args
    output["uri"] = uri

    return output
