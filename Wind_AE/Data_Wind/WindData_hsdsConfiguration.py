#Configure .hscfg to get wind from NREL
#verified it is work on 2023- Sept
import os

endpoint='https://developer.nrel.gov/api/hsds'
username=None
password=None
#api_key='3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf' #public
api_key='6DyGbHKzEhkNRVeC5XXRJ2WcBSBWrusoeis1LzG5'  #Victor

def saveConfig(username, password, endpoint, api_key):

    filepath = os.path.expanduser('~/.hscfg')
    print("Saving config file to: {}".format(filepath))
    with open(filepath, 'w') as file:
        file.write("# HDFCloud configuration file\n")
        if endpoint:
            file.write("hs_endpoint = {}\n".format(endpoint))
        else:
            file.write("hs_endpoint = \n")
        if username:
            file.write("hs_username = {}\n".format(username))
        else:
            file.write("hs_username = \n")
        if password:
            file.write("hs_password = {}\n".format(password))
        else:
            file.write("hs_password = \n")
        if api_key:
            file.write("hs_api_key = {}\n".format(api_key))
        else:
            file.write("hs_api_key = \n")
            
saveConfig(username, password, endpoint, api_key)
