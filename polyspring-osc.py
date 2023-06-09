from argparse import ArgumentParser
from polyspring import Corpus
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client
from shapely import Polygon
from math import ceil

class CorpusMax(Corpus):
    def __init__(self, track, cols, client):
        Corpus.__init__(self, track, cols)
        self.client = client

    def export(self, interp=0):
        current_idx = 0
        for key, length in self.buffers_md.items():
            self.client.send_message('/buffer_index', int(key))
            self.client.send_message('/matrixcol', 2)
            buffer = self.points[current_idx : current_idx + length]
            current_idx += length            
            uniX = [p.scaled_x * (1 - interp) + p.scaled_og_x * interp for p in buffer]
            uniY = [p.scaled_y * (1 - interp) + p.scaled_og_y * interp for p in buffer]
            n_rows = len(uniX)
            steps = int(ceil(n_rows/200))
            for i in range(steps):
                if i != steps-1:
                    self.client.send_message('/set_matrix', [i*200] + uniX[i*200:(i+1)*200])
                else :
                    self.client.send_message('/set_matrix', [i*200] + uniX[i*200:])
            self.client.send_message('/matrixcol', 3)
            for i in range(steps):
                if i != steps-1:
                    self.client.send_message('/set_matrix', [i*200] + uniY[i*200:(i+1)*200])
                else :
                    self.client.send_message('/set_matrix', [i*200] + uniY[i*200:])
        self.client.send_message('/refresh', 1)


# Functions mapped to OSC addresses
# ---- Manage import from Max
def import_init(addrs, args, *message):
    print('--> Export from Max...')
    args[1]['buffer'] = {}
    args[1]['osc_batch_size'] = int(message[1])
    args[1]['nb_buffer'] = int(message[0])
    args[1]['nb_lines'] = {}
    args[1]['remaining_lines'] = {}
    args[1]['cols'] = tuple()
    args[0].send_message('/begin_import', 1)

def add_buffer(addrs, args, *message):
    n_cols = int(message[2]) + 1
    n_rows = int(message[0])
    buffer = str(message[1])
    args[1]['buffer'][buffer] = [[0. for i in range(n_cols)]for j in range(n_rows)]
    args[1]['remaining_lines'][buffer] = n_rows
    args[1]['nb_lines'][buffer] = 0
    args[0].send_message('/start_dump', 1)

def add_line(addrs, args, *message):
    index = int(message[-2])
    buffer = str(message[-1])
    n_descr = len(message) - 2
    descriptors = [message[i] for i in range(n_descr)]
    args[1]['buffer'][buffer][index] = descriptors
    # update line count
    args[1]['remaining_lines'][buffer] -= 1
    args[1]['nb_lines'][buffer] += 1
    # check if all the buffer has been imported
    if args[1]['remaining_lines'][buffer] == 0:
        print('buffer', buffer, ',',args[1]['nb_lines'][buffer], 'grains' )
        args[0].send_message('/next_buffer', 1)
        end_test = [item<=0 for i,item in args[1]['remaining_lines'].items()]
        len_test = args[1]['nb_buffer'] == len(end_test)
        if all(end_test) and len_test:
            args[0].send_message('/done_import', 1)
    # ask for next batch to be sent if this one is done
    if args[1]['nb_lines'][buffer] % args[1]['osc_batch_size'] == 0:
        args[0].send_message('/next_batch', 1)

def set_cols(addrs, args, *cols):
    args[1]['corpus'].setCols(cols)

def write_track(addrs, args, *cols):
    xcol, ycol = cols
    for idx_buffer, track in args[1]['buffer'].items():
        args[0].send_message('/buffer_index', int(idx_buffer))
        for row in track:
            grain = [row[0], row[xcol], row[ycol], row[xcol], row[ycol]]
            args[0].send_message('/append', grain)
    args[0].send_message('/done_init', 1) 
    args[0].send_message('/update', 1)
    args[1]['corpus'] = CorpusMax(args[1]['buffer'], (xcol, ycol), args[0])
    print('<-- Done')


# ---- Manage polyspring
def distribute(addrs, args, *unused):
    print('--> Distributing...')
    c1, c2 = args[1]['corpus'].distribute(exportPeriod=1)
    args[1]['corpus'].export()
    if c1 < 0:
        print('<-- Force Stop ({} steps, {} triangulations)'.format(-c1, c2))
    else:
        print('<-- Done ({} steps, {} triangulations)'.format(c1, c2))
    client.send_message('/update', 1)

def change_interp(addrs, args, interp_value):
    args[1]['corpus'].export(float(interp_value))

def change_region(addrs, args, *coord):
    if args[1]['available']:
        print('--- change region')
        vertices = [(coord[i],1-coord[i+1]) for i in range(0,len(coord),2)]
        region = Polygon(vertices)
        args[1]["corpus"].setRegion(region)

def change_density(addrs, args, func):
    print('--- change density')
    args[1]["corpus"].h_dist = eval('lambda x, y :' + str(func))

def stop(addrs, args, *unused):
    args[1]["corpus"].stop_distribute()


# ---- Attractors
def attractors(addrs, args, *param):
    if len(param) % 5 == 0:
        gaussians_param = [param[5*i:5*(i+1)] for i in range(len(param)//5)]
        args[1]["corpus"].simple_attractors(gaussians_param)
    else:
        args[1]["corpus"].simple_attractors(gaussians_param, reset=True)


if __name__ == "__main__":
    print('Starting server...')
    # Client side (send to Max)
    parser_client = ArgumentParser()
    parser_client.add_argument("--ip", default="127.0.0.1")
    parser_client.add_argument("--port", type=int, default=8012)
    args_client = parser_client.parse_args()
    client = udp_client.SimpleUDPClient(args_client.ip, args_client.port)

    # Server side parameters (receive from Max)
    parser_server = ArgumentParser()
    parser_server.add_argument("--ip", default="127.0.0.1")
    parser_server.add_argument("--port", type=int, default=8011)
    args_server = parser_server.parse_args()

    # Init the global hash table and the dispatcher
    global_hash = {'buffer':{}, 'available':False}
    dispatcher = dispatcher.Dispatcher()

    # Map OSC adresses to functions
    # Manage import from Max ----
    dispatcher.map("/export_init", import_init, client, global_hash)
    dispatcher.map("/add_buffer", add_buffer, client, global_hash)
    dispatcher.map("/add_line", add_line, client, global_hash)
    dispatcher.map("/set_cols", set_cols, client, global_hash)
    dispatcher.map("/write_track", write_track, client, global_hash)
    # Manage unispring ----
    dispatcher.map("/distribute", distribute, client, global_hash)
    dispatcher.map("/interpolation", change_interp, client, global_hash)
    dispatcher.map("/region", change_region, client, global_hash)
    dispatcher.map("/density", change_density, client, global_hash)
    dispatcher.map("/attractors", attractors, client, global_hash)
    dispatcher.map("/stop", stop, client, global_hash)

    # Init server
    server = osc_server.ThreadingOSCUDPServer((args_server.ip, args_server.port), dispatcher)

    # Launch the server
    print("----- Serving on {}".format(server.server_address))
    server.serve_forever()