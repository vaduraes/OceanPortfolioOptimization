import numpy as np

def compute_distance(latlong1, latlong2):
    """
    Compute the distance between a latlong value and a vector of latlong values.
    
    Arguments:
    latlong1 -- Tuple representing the latlong value.
    latlong2 -- List of tuples representing the latlong values of the vector.
    
    Returns:
    List of distances between the latlong value and each point in the vector.
    """
    lat1, lon1 = latlong1
    lat2, lon2 = np.array(latlong2).T
    
    # Convert degrees to radians
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)
    
    # Haversine formula
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad
    a = np.sin(delta_lat/2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(delta_lon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    # Radius of the Earth in kilometers
    r = 6371
    
    # Calculate the distance
    distances = r * c
    
    return distances

#latlong1: Transmission
#latlong2: Devices
#latlong1 is a vector of points and latlong2 is a vector of points
def GetIdxOutRadious(latlong1, latlong2, Radious):
    IdxOut=[]
    #if latlong2=[], return IdxOut=[]
    try:
        for i in range(latlong1.shape[0]):
            D=compute_distance(latlong1[i,:], latlong2)
            IdxOut.append(np.where(D>=Radious)[0])
        return IdxOut
    except:
        if len(latlong2)==0:
            IdxOut=[[]]*latlong1.shape[0]
            return IdxOut
        else:
            return print("Error in GetIdxOutRadious function")


#Which turbines are inside the radious of the tranmission system
#latlong1 is a vector of points and latlong2 is a vector of points
def GetIdxInRadious(latlong1, latlong2, Radious):
    IdxOut=[]
    
    try:
        for i in range(latlong1.shape[0]):
            D=compute_distance(latlong1[i,:], latlong2)
            IdxOut.append(np.where(D<=Radious)[0])
        return IdxOut
    except:
        if len(latlong2)==0:
            IdxOut=[[]]*latlong1.shape[0]
            return IdxOut
        else:
            return print("Error in GetIdxInRadious function")


#In this case we determine the index of the devices that are inside the radious of the transmission system
#Considering latlong1 as a single point and latlong2 as a vector of points
def GetIdxInRadious_Simple(latlong1, latlong2, Radious):
    IdxIn=[]
    D=compute_distance(latlong1, latlong2)
    return np.where(D<=Radious)[0]