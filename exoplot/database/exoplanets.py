import sqlite3

def load_payload_from_database(ra_start,ra_end,dec_start,dec_end, args):
    payload = []
    # todo: load the database from within the distribution
    # alternatively... use --database argument and load it from outside the package (makes it easier to update)
    # https://stackoverflow.com/questions/25270817/loading-sqlite-database-in-python-distribution
    dbname = "exoplanets.sqlite3"

    conn = sqlite3.connect(dbname)
    query = f"SELECT ra, dec, hostname, sy_pnum FROM exoplanets WHERE ra>{ra_start}-5 AND ra<{ra_end}+5 " \
            f"AND dec>{dec_start}-5 AND dec<{dec_end}+5"

    # very expensive to try to plot all the exoplanets, but it saves the complexity of doing a real cone search
    #query = f"SELECT ra, dec, hostname, sy_pnum FROM exoplanets"
    print(query)

    cursor = conn.execute(query)
    size1 = int(args.size)
    size2 = size1 + 10
    for row in cursor:
        record = {"ra": row[0], "dec": row[1], "label": row[2], "shape": "exoplanet", "size": size1, "color": "red"}
        payload.append(record)

        # if this star as multiple exoplanets, then also draw a green circle
        if row[3]> 1:
            record = {"ra": row[0], "dec": row[1], "": row[2], "shape": "exoplanet", "size": size2, "color": "green"}
            payload.append(record)

    print(f"{len(payload)} exoplanets found")

    conn.close()

    return payload