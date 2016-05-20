import iris

class constants():

    def __init__(self):

        self.Mw=iris.cube.Cube(18.0160, units='g mol**-1',
                    long_name='molecular weight of water')
        self.Md=iris.cube.Cube(28.9660, units='g mol**-1',
                    long_name='molecular weight of dry air')
        self.R = iris.cube.Cube(8.31432E3, units='J mol**-1 K**-1',
                    long_name='gas constant')
        self.cpd  = iris.cube.Cube(1005.7, units='J kg**-1 K**-1',
                    long_name='specific heat of dry air at constant pressure')
        self.cpv = iris.cube.Cube(1864., units='J kg**-1 K**-1',
                    long_name='specific heat of water vapour at constant pressure and 300K')
        self.rhod = iris.cube.Cube(1.2, units='kg m**-3',
                    long_name='specific mass of dry air for standard atmosphere')
        self.Rd = iris.cube.Cube(self.R.data/self.Md.data, units='J kg**-1 K**-1',
                    long_name='specific gas constant for dry air')
        self.Rv = iris.cube.Cube(self.R.data/self.Mw.data, units='J kg**-1 K**-1',
                    long_name='specific gas constant for water vapour')
        self.eps = self.Mw/self.Md
        self.eps.long_name= 'ratio of the molecular weights of water and dry air'
        self.g = iris.cube.Cube(9.81, units='m s**-2',
                    long_name='gravitational acceleration')


