import argparse
import unispring_v2 as usp
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client
import shapely as sh
import numpy as np

def MinMaxScale(track):
    n_descr = len(track['1'][0])
    norm_track = {}
    list_min = [float('inf') for i in range(n_descr)]
    list_max = [float('-inf') for i in range(n_descr)]
    for key, table in track.items():
        for line in table:
            for i in range(1,n_descr):
                if line[i] < list_min[i]:
                    list_min[i] = line[i]
                elif line[i] > list_max[i]:
                    list_max[i] = line[i]
    for key, table in track.items():
        norm_track[key] = []
        for line in table:
            new_line = [line[0]]
            for i in range(1,n_descr):
                new_line.append((line[i]-list_min[i])/(list_max[i]-list_min[i]))
            norm_track[key].append(new_line)
    return norm_track

# Functions mapped to OSC addresses
# ---- Manage import from Max
def import_init(addrs, args, *message):
    print('--> Export from Max...')
    args[1]['buffer'] = {}
    args[1]['osc_batch_size'] = int(message[1])
    args[1]['nb_buffer'] = int(message[0])
    args[1]['nb_lines'] = {}
    args[1]['remaining_lines'] = {}
    args[0].send_message('/begin_import', 1)

def add_buffer(addrs, args, *message):
    n_cols = int(message[2]) + 1
    n_rows = int(message[0])
    buffer = str(message[1])
    args[1]['buffer'][buffer] = np.zeros((n_rows, n_cols))
    args[1]['remaining_lines'][buffer] = n_rows
    args[1]['nb_lines'][buffer] = 0
    args[0].send_message('/start_dump', 1)

def add_line(addrs, args, *message):
    index = int(message[-2])
    buffer = str(message[-1])
    n_descr = len(message) - 2
    descriptors = np.asarray([message[i] for i in range(n_descr)])
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

def create_norm_track(addrs, args, *unused):
    args[1]['norm_buffer'] = {}
    args[1]['norm_buffer'] = MinMaxScale(args[1]['buffer'])
    args[0].send_message('/done_create', 1)

def write_norm_track(addrs, args, *unused):
    for idx_buffer, track in args[1]['norm_buffer'].items():
        args[0].send_message('/buffer_index', int(idx_buffer))
        for i,line in enumerate(track):
            args[0].send_message('/append', line)
    args[0].send_message('/done_norm', 1) 
    args[0].send_message('/update', 'update')
    args[1]['corpus'] = usp.Corpus(args[1]['norm_buffer'], args[0])
    args[1]['available'] = True
    print('<-- Done')

# ---- Manage unispring
def unispring(addrs, args, *descr):
    if args[1]['available']:
        print('--> Distributing...')
        args[1]['available'] = False
        args[1]['corpus'].setDescr(descr[0]+1, descr[1]+1)
        c1, c2 = args[1]['corpus'].unispring(exportPeriod=1, switch_on=args[1]['available'])
        args[1]['corpus'].exportToMax()
        args[1]['available'] = True
        print('<-- Done ({} steps, {} triangulations)'.format(c1, c2))

def change_interp(addrs, args, interp_value):
    args[1]['corpus'].setInterp(float(interp_value))
    args[1]['corpus'].exportToMax()

def change_region(addrs, args, *coord):
    if args[1]['available']:
        print('--- change region')
        vertices = [(coord[i],1-coord[i+1]) for i in range(0,len(coord),2)]
        region = sh.Polygon(vertices)
        args[1]["corpus"].setRegion(region)

def change_density(addrs, args, func):
    print('--- change density')
    args[1]["corpus"].h_dist = eval('lambda x, y :' + str(func))


# ---- Attractors
def attractors(addrs, args, *param):
    if args[1]['available']:
        args[1]['available'] = False
        if len(param) % 5 == 0:
            gaussians_param = [param[5*i:5*(i+1)] for i in range(len(param)//5)]
            args[1]["corpus"].simple_attractors(gaussians_param)
        else:
            args[1]["corpus"].simple_attractors(gaussians_param, reset=True)
        args[1]['available'] = True


if __name__ == "__main__":
    print('Starting server...')
    # Client side (send to Max)
    parser_client = argparse.ArgumentParser()
    parser_client.add_argument("--ip", default="127.0.0.1")
    parser_client.add_argument("--port", type=int, default=8012)
    args_client = parser_client.parse_args()
    client = udp_client.SimpleUDPClient(args_client.ip, args_client.port)

    # Server side parameters (receive from Max)
    parser_server = argparse.ArgumentParser()
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
    dispatcher.map("/create_norm_track", create_norm_track, client, global_hash)
    dispatcher.map("/write_norm_track", write_norm_track, client, global_hash)
    # Manage unispring ----
    dispatcher.map("/unispring", unispring, client, global_hash)
    dispatcher.map("/interpolation", change_interp, client, global_hash)
    dispatcher.map("/region", change_region, client, global_hash)
    dispatcher.map("/density", change_density, client, global_hash)
    dispatcher.map("/attractors", attractors, client, global_hash)

    # Init server
    server = osc_server.ThreadingOSCUDPServer((args_server.ip, args_server.port), dispatcher)

    # Launch the server
    print("----- Serving on {}".format(server.server_address))
    server.serve_forever()