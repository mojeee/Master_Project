if __name__ == '__main__':
    # import required libraries
    import h5py as h5
    import numpy as np
    import matplotlib.pyplot as plt
    

# open the file
    f=h5.File("prof000009.h5",'r')
  # get the keys
    #key_list = f.keys()
    #print key_list
    #with f['T'].value we can print data
    key_list='T'
    data=f[key_list]
    dataf=data[:]
    key_list='x'
    data=f[key_list]
    dataf2=data[:]
    f.close()
    plt.plot(dataf2, dataf,'r')
    plt.title("Laminar Premix Flame")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Position(m)')
    
    plt.show()
    plt.savefig('test.png')
