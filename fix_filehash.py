import os.path
import getpass
from flows.aadc_db import AADC_DB
from flows.utilities import get_filehash
from tendrils.utils import load_config


if __name__ == '__main__':
    config = load_config()

    with AADC_DB() as db:

        db.cursor.execute("SELECT * FROM flows.files WHERE targetid=8;")

        for row in db.cursor.fetchall():
            if row['path'].endswith('.fits.gz'):
                continue

            p = row['path'].replace('.fits', '.fits.gz')

            fpath = os.path.join('/data/flows/archive/archive', p)
            try:
                fsize = os.path.getsize(fpath)
                fhash = get_filehash(fpath)
                print(p)
                print(fsize)
                print(fhash)

                db.cursor.execute("UPDATE flows.files SET path=%s,filehash=%s,filesize=%s WHERE fileid=%s;", [p, fhash, fsize, row['fileid']])
                db.conn.commit()
                #break
            except FileNotFoundError:
                print(f'File {fpath} not found')
