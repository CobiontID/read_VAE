#!/usr/bin/env python

"""
<cw21@sanger.ac.uk>

A minimal server that receives a fasta file and passes it to blast via a subprocess.
When the subprocess returns, it sends back the results in the format specified in the variable 'out_fmt'.

Usage: To launch the server, run 'python minimal_blast_server.py'. Set PORT to an open port.
From a remote machine, run 'curl -T query.fa http://<IP>:<PORT>' where IP is the IP address of this machine.

Loading the database into shared memory with 'mkdir -p /dev/shm/nt; parallel cp {} /dev/shm/nt/ ::: /path_to_db/nt_2021_06/*'
ensures faster turnaround for each individual query.
"""

import logging
import tempfile
import os

try:
    from urllib.parse import unquote
except ImportError:
    # Python 2.
    from urllib import unquote

import tornado.ioloop
import tornado.web
from tornado import options
import subprocess


# Configuration
#
# Port to listen on
PORT = 35227
# Number of threads for blast
NTHREADS = 8

# Path to blast executable
blastn_path = "/usr/bin/blastn"
# Path to blast database
nt_db_path = "/dev/shm/nt/nt"
# Output format
out_fmt = "6 std staxid stitle"



def format_query(file_path):
    # TODO: Check query is formatted correctly first?
    return f'{blastn_path} -db {nt_db_path} -query {file_path} -num_threads {NTHREADS} -outfmt "{out_fmt}"'

def blast(file_path):
    query_string = format_query(file_path)
    blast = subprocess.run(query_string, capture_output=True, shell=True)
    if blast.returncode == 0:
        if len(blast.stdout) == 0:
             return "No hits returned\n"
        else:
              # Filtering of results could be added here
              return "{}".format(blast.stdout.decode('utf-8'))
    else:
        return "Error, {}".format(blast)

class POSTHandler(tornado.web.RequestHandler):
    def post(self):
        #temp = open("upload.file", "wb")
        for field_name, files in self.request.files.items():
            for info in files:
                filename, content_type = info["filename"], info["content_type"]
                body = info["body"]
                #temp.write(body)
                logging.info(f'POST "{filename}" "{content_type}" {len(body)} bytes')
                self.write("OK")

        self.write("OK")


@tornado.web.stream_request_body
class PUTHandler(tornado.web.RequestHandler):
    
    def initialize(self):
        self.temp = None
        self.bytes_read = 0

    def data_received(self, chunk):
        self.bytes_read += len(chunk)
        if self.temp is None:
            self.temp = tempfile.NamedTemporaryFile(delete=False)
        self.temp.write(chunk)

    def put(self, filename):
        filename = unquote(filename)
        mtype = self.request.headers.get("Content-Type")
        logging.info(f'PUT "{filename}" "{mtype}" {self.bytes_read} bytes')
        
        if self.temp is not None:
            self.temp.flush()
            os.fsync(self.temp.fileno())
            print(self.temp.name)
            matches = blast(self.temp.name)

            self.temp.close()
            os.unlink(self.temp.name)
            self.temp = None

            self.write(matches)

def make_app():
    return tornado.web.Application([(r"/post", POSTHandler), (r"/(.*)", PUTHandler)])


if __name__ == "__main__":
    # Tornado configures logging.
    options.parse_command_line()
    app = make_app()
    app.listen(PORT)
    tornado.ioloop.IOLoop.current().start()

