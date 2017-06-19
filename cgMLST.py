import os
import sys
import shutil
import argparse


class ENAUploader():

    def __init__(self):

        self.ftp_upload = None


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="")
    # Posotional arguments
    parser.add_argument("input",
                        help="Tab seperated text file containing all the\
                              metadata needed.",
                        metavar='INPUT',
                        default=None)
    # Optional arguments
    parser.add_argument("-o", "--output",
                        help="Output prefix.",
                        default="./metadata",
                        metavar="OUTPUT_PREFIX")
    parser.add_argument("-c", "--credentials",
                        help="File with user name and password for ENA\
                              account. If left out user must input credentials\
                              via the --ena_user and --ena_pass flags.",
                        default=None,
                        metavar="TXT_FILE")
    parser.add_argument("-eu", "--ena_user",
                        help="ENA user name. Not necessary if credentials file\
                              has been provided.",
                        default=None,
                        metavar="STRING")
    parser.add_argument("-ep", "--ena_pass",
                        help="ENA password. Not necessary if credentials file\
                              has been provided.",
                        default=None,
                        metavar="STRING")
    parser.add_argument("--center",
                        help="Sequencing center. Default is DTU-GE.",
                        default="DTU-GE",
                        metavar="STRING")
    parser.add_argument("--title",
                        help="Title of project/study.",
                        default=None,
                        metavar="STRING")
    parser.add_argument("--abstract",
                        help="Abstract/Description of project/study.",
                        default=None,
                        metavar="STRING")
    parser.add_argument("--material",
                        help="TODO: Find description on ENA website.",
                        default="DNA",
                        metavar="STRING")
    parser.add_argument("--selection",
                        help="TODO: Find description on ENA website.",
                        default="genome",
                        metavar="STRING")
    parser.add_argument("--broker",
                        help="Name of broker account. Default is DTU-GE. If no\
                              broker is needed, an empty string should be\
                              provided.",
                        default="DTU-GE",
                        metavar="STRING")
    parser.add_argument("--umbrella",
                        help="Umbrella accession number. Default is PRJEB6071\
                              (Acc. number used for all CGE uploads to ENA).\
                              If no umbrella accession number is needed, an\
                              empty string should be provided.",
                        default="PRJEB6071",
                        metavar="STRING")
    parser.add_argument("--type",
                        help="Specify if the uploaded data is isolate data or\
                              metagenomic data.",
                        choices=["isolate", "metagenomic"],
                        default=None)
    parser.add_argument("--release",
                        help="Specify release date. Format: YYYY-MM-DD.",
                        default=None,
                        metavar="DATE")
    parser.add_argument("--seperator",
                        help="The input table is per default assumed to be\
                              seperated by tabs. This option makes it\
                              possible to use other operators like commas,\
                              semicolons etc. However, a seperator must be\
                              a single character.",
                        default="\t",
                        metavar="STRING")
    parser.add_argument("--upload_method",
                        help="Specify which method you want to use to upload\
                              sequence data to ENA. Default is FTP. IMPORTANT:\
                              ASPERA has not yet been implemented.",
                        choices=["FTP", "ASPERA"],
                        default="FTP")
    parser.add_argument("--upload_limit",
                        help="This option is only by the ASPERA upload method.\
                              Set maximum bandwith aspera is allowed to use.\
                              Default is 20 MB/s.",
                        metavar="MB/s",
                        type=int,
                        default=20)
    parser.add_argument("--checklist",
                        help="Specify which ENA checklist to test against.\
                              Default is 'gmi'",
                        choices=["gmi", "minimum"],
                        default="gmi")
    parser.add_argument("--test",
                        help="Submission to test server.",
                        action="store_true",
                        default=False)

    args = parser.parse_args()

    # Check credentials file.
    if args.credentials:
        ena_auth = Auth(args.credentials)
    elif args.ena_user and args.ena_pass:
        ena_auth = Auth()
        ena_auth.user = args.ena_user
        ena_auth.password = args.ena_pass
    else:
        print("ERROR: Credentials file not found and no ENA user name or\
               password was provided.")
        quit(1)

    if not args.title:
        print("ERROR: Title can not be left empty.")
        quit(1)

    if not args.abstract:
        print("ERROR: Abstract can not be left empty.")
        quit(1)

    if not args.type:
        print("ERROR: Type can not be left empty.")
        quit(1)

    if not args.release:
        print("ERROR: Release date can not be left empty.")
        quit(1)

    if not args.broker:
        args.broker = None

    if not args.umbrella:
        args.umbrella = None

    #
    # Extract data from input and create objects
    #

    prj_data = {
        "ena_user": ena_auth.user,
        "ena_pass": ena_auth.password,
        "center": args.center,
        "prj_title": args.title,
        "prj_abstract": args.abstract,
        "prj_type": args.type,
        "release": args.release,
        "material": args.material,
        "selection": args.selection,
        "umbrella": args.umbrella,
        "broker": args.broker
    }

    ena_project = ENAUploader(**prj_data)

    input_data = []

    with open(args.input, "r") as fh:
        header_line = fh.readline()
        header_line = header_line.rstrip()
        headers = header_line.split(args.seperator)

        header_column = {}
        for i, header in enumerate(headers):
            header_column[header] = i

        for line in fh:
            line_data = {}

            line.rstrip()
            if not line:
                continue

            line_list = line.split(args.seperator)
            # Append None to the end of the list. If a header is not
            # found in the input data it will retrieve the last
            # element (i.e. None).
            line_list.append(None)

            #
            # Extract data from each recognized column
            #
            # Note: Column header names are defined as class static class
            # variables in the ENAUploader class

            line_data[ENAUploader.voc_seq_plat] = line_list[
                header_column.get(ENAUploader.voc_seq_plat, -1)]

            table_layout = line_list[
                header_column.get(ENAUploader.voc_seq_type, -1)]

            if table_layout:
                table_layout = table_layout.upper()

            if table_layout != "PAIRED" and table_layout != "SINGLE":
                print("ERROR: Sequecning layout/type must be either 'SINGLE' "
                      "or 'PAIRED'")
                quit(1)

            line_data[ENAUploader.voc_seq_type] = table_layout

            line_data[ENAUploader.voc_instrument] = line_list[
                header_column.get(ENAUploader.voc_instrument, -1)]

            line_data[ENAUploader.voc_insert] = line_list[
                header_column.get(ENAUploader.voc_insert, -1)]

            line_data[ENAUploader.voc_file1] = line_list[
                header_column.get(ENAUploader.voc_file1, -1)]

            line_data[ENAUploader.voc_file2] = line_list[
                header_column.get(ENAUploader.voc_file2, -1)]

            line_data[ENAUploader.voc_uuid] = line_list[
                header_column.get(ENAUploader.voc_uuid, -1)]

            line_data[ENAUploader.voc_sample_name] = line_list[
                header_column.get(ENAUploader.voc_sample_name, -1)]

            line_data[ENAUploader.voc_isolate] = line_list[
                header_column.get(ENAUploader.voc_isolate, -1)]

            line_data[ENAUploader.voc_strain] = line_list[
                header_column.get(ENAUploader.voc_strain, -1)]

            line_data[ENAUploader.voc_org] = line_list[
                header_column.get(ENAUploader.voc_org, -1)]

            table_lat = line_list[
                header_column.get(ENAUploader.voc_lat, -1)]
            if table_lat:
                table_lon = line_list[
                    header_column.get(ENAUploader.voc_lon, -1)]
            else:
                table_lon = None

            line_data[ENAUploader.voc_lat] = table_lat
            line_data[ENAUploader.voc_lon] = table_lon

            table_country = line_list[
                header_column.get(ENAUploader.voc_country, -1)]

            line_data[ENAUploader.voc_country] = table_country

            table_region = line_list[
                header_column.get(ENAUploader.voc_region, -1)]

            line_data[ENAUploader.voc_region] = table_region

            table_city = line_list[
                header_column.get(ENAUploader.voc_city, -1)]

            line_data[ENAUploader.voc_city] = table_city

            line_data[ENAUploader.voc_col_date] = line_list[
                header_column.get(ENAUploader.voc_col_date, -1)]

            table_host = line_list[
                header_column.get(ENAUploader.voc_host, -1)]
            table_host_stat = line_list[
                header_column.get(ENAUploader.voc_host_stat, -1)]
            if table_host and table_host_stat:
                table_host_stat = table_host_stat.lower()
                if (table_host_stat and
                    table_host_stat != "healthy" and
                    table_host_stat != "diseased" and
                    table_host_stat != "not applicable" and
                    table_host_stat != "not collected" and
                    table_host_stat != "not provided" and
                        table_host_stat != "restricted access"):
                    print("ERROR: If a host is provided. The status of the"
                          "host must either be 'healthy', 'diseased', 'not "
                          "applicable', 'not collected', 'not provided', or "
                          "'restricted access'")
                    quit(1)

            line_data[ENAUploader.voc_host] = table_host
            line_data[ENAUploader.voc_host_stat] = table_host_stat

            line_data[ENAUploader.voc_col_by] = line_list[
                header_column.get(ENAUploader.voc_col_by, -1)]

            line_data[ENAUploader.voc_iso_src] = line_list[
                header_column.get(ENAUploader.voc_iso_src, -1)]

            line_data[ENAUploader.voc_serovar] = line_list[
                header_column.get(ENAUploader.voc_serovar, -1)]

            input_data.append(line_data)

    try:
        ena_project.data(input_data)
    # Object creation requires an internet connection. A missing connection
    # will raise a requests.exceptions.ConnectionError.
    except requests.exceptions.ConnectionError as err:
        print("Submission failed.")
        print("Please check your internet connection.")
        quit(1)

    # Check input data and report errors
    report = ena_project.validate(checklist=args.checklist)
    # The number of errors are recorded in the returned object.
    if report.error_count:
        print(report.report_to_str())
        print("!! Metadata contained errors.")
        print("!! Submission has been cancelled.")
        quit(1)

    print("Validation successful.")
    print("Uploading...")

    if args.upload_method == "ASPERA":
        print("ERROR: ASPERA has not yet been implemented.")
        print("Switching to FTP")
        args.upload_method = "FTP"

    if args.upload_method == "FTP":
        ena_project.upload("FTP")
        print("Uploading... done!")

    print("Submitting to ENA...")

    try:
        ena_response = ena_project.submit_xmls(output=args.output,
                                               test_server=args.test)
    except OSError as err:
        print(err.args[0])
        uploader_obj.clean()
        quit(1)

    #
    # Handle ENA response & print accessions
    #

    if ena_response.success:
        print("Submission successful.")

        submission_accession = ena_response.accessions["SUBMISSION"]
        print("Submission accession no.: " + submission_accession)

        project_accession = ena_response.accessions["PROJECT"]
        print("Project accession no.: " + project_accession)

        accession_output = "Submission acc.: " + submission_accession + "\n"
        accession_output += "Project acc.: " + project_accession + "\n\n"

        accession_output += "ID\tSample\tExperiment\tRun\n"

        for line_data in input_data:
            # Get sample information from your own database
            sample_uuid = line_data[ENAUploader.voc_uuid]

            # Retrieve the relevant accession numbers for the sample
            accs = ena_project.get_rel_sample_acc(
                uuid=sample_uuid, ena_response=ena_response)

            # There is only one sample accession no. per sample
            accession_output += sample_uuid + "\t" + accs["SAMPLE"] + "\t"

            # There can be multiple experiment accessions per sample
            # Most of the time there will be just one
            accession_output += ",".join(accs["EXPERIMENT"]) + "\t"

            # There can be multiple run accessions per sample
            # Most of the time there will be just one
            accession_output += ",".join(accs["RUN"]) + "\n"

        with open(args.output + "_acc_numbers.txt", "w") as fh:
            fh.write(accession_output)

        print("Wrote: " + args.output + "_acc_numbers.txt")
    else:
        for message in ena_response.messages:
            print(message)
        print("!! Submission failed.")
