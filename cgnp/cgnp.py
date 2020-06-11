import mbuild as mb
import numpy as np


class cgnp(mb.Compound):
    
    '''
    Builds a coarse-grained, silica, alkane-coated nanoparticle.

    Parameters
    ----------
    radius : float
        Radius of the nanopatile in nanometers
    bead_diameter : float
        Diameter fo the coarse-grained nanoparticle core in nanometers
    chain_length : mb.Compound
        Number of alkanes in single chain attached to nanoparticle surface. Must be a factor of three (cg ratio 3:1)
    chain_density : float
        Density of coating on nanoparticle surface in chains/nanometer squared
    backfill : mb.Compound, optional, default=None
        Place chains on uncoated sections of nanoparticle surface
    coating_pattern : str, optional, default='isotropic'
        TODO: Add support for coating patterns
    '''

    def __init__(self, radius=10, bead_diameter=0.2, chain_length=12, chain_density=0.5, backfill=None):
        super(cgnp, self).__init__()
       
        
        ''' Makes the core '''
        class Core(mb.Compound):
            def __init__ (self, radius):
                super(Core, self).__init__()
                
                core_bead = mb.Particle(name='_CGN')

                # calculate n number of core particles from radius and surface ratio
                n_core_particles = (int)(0.2 * (4*np.pi*(radius**2)))
                
                # arrange beads into a sphere pattern and add it to the compound
                core_pattern = mb.SpherePattern(n_core_particles)
                core_pattern.scale(radius)

                core_formation = core_pattern.apply(core_bead)
                self.add(core_formation)
                print('Core beads added.')

                # apply ports to each bead in core_formation
                for i, pos in enumerate(core_pattern.points):
                    port = mb.Port(anchor=core_formation[i], orientation=pos, separation=radius/5)
                    self.add(port, label='port[$]')

                print('Ports added to core')

        ''' Makes the coarse-grained alkane chain '''
        class CGAlkane(mb.Compound):
            def __init__ (self, chain_length):
                super(CGAlkane, self).__init__()
                
                if (chain_length%3 != 0):
                    raise Exception ('Due to 3:1 cg-ratio of alkane chains on this nanoparticle model, chain length must be divisible by three.')
                chain_length = (int)(chain_length/3)
                chain_separation = 0.30

                ''' Makes the beads to populate the middle of the chain '''
                class CGMMM(mb.Compound):
                    def __init__ (self, chain_separation=0.30):
                        super(CGMMM, self).__init__()
                        middle_bead = mb.Particle(name='_MMM')
                        self.add(middle_bead, 'middle_bead')
                        
                        self.add(mb.Port(anchor=self['middle_bead']), label='up')
                        self.add(mb.Port(anchor=self['middle_bead']), label='down')
                        mb.Compound.translate(self['up'], [0, chain_separation, 0])
                        mb.Compound.translate(self['down'], [0, -(chain_separation), 0])
                
                ''' Makes the beads to populate the end of the chain. Should only have one per alkane chain '''
                class CGMME(mb.Compound):
                    def __init__ (self, chain_separation=0.30):
                        super(CGMME, self).__init__()
                        end_bead = mb.Particle(name='_MME')
                        self.add(end_bead, 'end_bead')
                        
                        self.add(mb.Port(anchor=self[0]), label='end')
                        mb.Compound.translate(self['end'], [0, chain_separation, 0])

                last_bead = CGMMM()
                self.add(last_bead)
                for i in range (chain_length-2):
                    current_bead = CGMMM()
                    mb.force_overlap(move_this=current_bead, from_positions=current_bead['up'], to_positions=last_bead['down'])
                    self.add(current_bead)
                    last_bead = current_bead

                end_bead = CGMME()
                self.add(end_bead)
                mb.force_overlap(move_this=end_bead, from_positions=end_bead['end'], to_positions=current_bead['down'])
                print('Chain successfully built.')


        np_core = Core(radius)
        self.add(np_core)
        print('Core added to Compound')

        n_cgn = self.n_particles
        i = n_cgn
        
        while (i > 0):
            alkane = CGAlkane(chain_length)
            self.add(alkane)
            hanging_port = alkane.all_ports()
            port_pos = hanging_port[0]
            mb.force_overlap(move_this=alkane, from_positions=port_pos, to_positions=np_core['port'][i-1])
            i -= 1
            print('Chain added. %d left.' % i)

        # TODO: Add support for multiple coating patterns
        # TODO: Add backfill support
        
        print('Finished adding chains. Now adding ridgid body labels.')
        self.label_rigid_bodies(rigid_particles='_CGN')
        
