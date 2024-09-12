# physics_numerical_calculation
# Memc Instructions

The module contains several functions to perform the numerical calculations related to mechanics and electromagnetism. The functions can calculate the mass, momentum of inertia, charge distribution, electric potential, electric field, magnetic vector potential and magnetic field for arbitrary objects depicted by the hull and shape3d models in mod3d module.

## Installation

To use this module, ensure you have the following dependencies installed:

```bash
pip install numpy scipy
```

The main module file is memc.py. You may only import the file only to start using it.
```python
import memc
```

The module relies on 3 other module files of volume_integration.py, conductor_potential.py and magnetic_integration. Please make sure they are all in correctly installed.

## Functions

### 1. `mass(hulls_dict)`
Calculates the mass/charge of objects by integrating their volume.

- **Parameters**:
  - `hulls_dict`: A dictionary to describe the hull objects.
                    hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
  
- **Returns**: 
  - The masses of the objects in the same dictionary format.
    {'object_name':
        {'integration': integration result (mass or charge), # the mass/charge of the object.
        'error': calculation error} # the integration error.
    }

```python
mass(hulls_dict)
```

### 2. `momentum(hulls_dict)`
Calculates the tensors of moment of inertia of objects.

- **Parameters**:
  - `hulls_dict`: A dictionary to describe the hull objects.
                    hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass density
                                        },
                                        ...
                                    }

- **Returns**:
  - The dictionary of momentum of inertia tensors.
    {'object_name': array([[ Ixx,  Ixy, Ixz],
                           [ Iyx,  Iyy, Iyz],
                           [ Izx,  Izy, Izz]])
     ...
    }

```python
momentum(hulls_dict)
```

### 3. `dipole(hulls_dict, point, electric = True)`
Calculates the dipole vectors of charged objects to a point.

- **Parameters**:
  - `hulls_dict`: A dictionary to describe the hull objects.
                    hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g. sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
  - `point`: the dipole reference point, (x, y, z): tuple, list or 1darray.
  - `electric`: if True, it calculates the electric dipole, else it calculates the gravitational dipole.

- **Returns**:
  - The dictionary of dipoles of objects.
    {'object_name': array([dipole in the x-direction, dipole in the y-direction, dipole in the z-direction]),
     ...
    }

```python
diple(hulls_dict, (1, 2, 3))
```

### 4. `potential(hulls_dict, point, electric = True)`
Calculates the electric/gravitational potentials of charged objects at a point.

- **Parameters**:
  - `hulls_dict`: A dictionary to describe the hull objects.
                    hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
 - `point`: the destination point where you want to calculate the electric field, (x, y, z): tuple, list or 1darray.
 - `electric`: if True, it calculates the electric potential, else it calculates the gravitational potential.
- **Returns**:
  - The dictionary of the electric potentials of objects contributing to the point.
    {'object_name':
        {'integration': potential at the point, # the electric potential to the destination point.
         'error': calculation error # the integration error.
         },
        ...
    }

```python
potential(hulls_dict, (1, 2, 3))
```

### 5. `field(hulls_dict, point, electric = True)`
Calculates the electric fields of charged objects at a point.

- **Parameters**:
  - `hulls_dict`: A dictionary to describe the hull objects.
                    hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g. sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
 - `point`: the destination point where you want to calculate the electric field, (x, y, z): tuple, list or 1darray.
 - `electric`: if True, it calculates the electric field, else it calculates the gravitational field.
- **Returns**:
  - The dictionary of the electric fields of objects contributing to the point.
    {'object_name': array([Ex, Ey, Ez]), # electric field vector of the object
     ...
    }

```python
field(hulls_dict, (1, 2, 3))
```

### 6. `conductor_charge(hulls_dict, division = 0.1)`
Calculates the charge distribution on conductors at given potentials. It has taken into account of the influences of other conductors.

- **Parameters**:
  - `hulls_dict`: A dictionary to describe the hull objects.
                    hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g. sphere()
                                        'potential': e.g. 1 # electric potential of the object
                                        },
                                        ...
                                    }
 - `division`: the spacing between adjacent points of calculation, default = 0.1.
- **Returns**:
  - The dictionary of the charge distributions of objects.
    {'object_name':
                    {'charge': 1darray of the charge on every point,
                    'points': {array[x, y, z], ...} # the coordinates (x, y, z) of every point of calculation on the surface of the conductor.
                    'total': total charge on the conductor
                    },
     ...
    }

```python
conductor_charge(hulls_dict)
```

### 7. `conductor_potential(charge_dist, point, epsilon = 0.1)`
Calculates the electric potentials contributed by conductors from the results of conductor_charge().

- **Parameters**:
  - `charge_dist`: the charge distribution of a conductor resulted from charge_distribution() calculation
  - `point`:  the destination point where you want to calculate the electric potential, (x, y, z): tuple, list or 1darray.
  - `epsilon`: a small value to prevent singularity of very small distance, default is 0.1. It is recommended to be equal to the value of the division of the hulls.
- **Returns**:
  - The tuple of the potential contributions
    ( { 'object_name': potential,
                ...
       }, # potential contributed by each object

       total_potential
    )

```python
conductor_potential(charge_dist, (1, 2, 3))
```

### 8. `conductor_field(charge_dist, point, epsilon = 0.1)`
Calculates the electric fields contributed by conductors from the results of conductor_charge().

- **Parameters**:
  - `charge_dist`: the charge distribution of a conductor resulted from charge_distribution() calculation
  - `point`:  the destination point where you want to calculate the electric field, (x, y, z): tuple, list or 1darray.
  - `epsilon`: a small value to prevent singularity of very small distance, default is 0.1. It is recommended to be equal to the value of the division of the hulls.
- **Returns**:
  - The tuple of the electric field contributions
    ( { 'object_name': field,
                ...
       }, # electric field contributed by each object

       total_electric_field
    )

```python
conductor_field(charge_dist, (1, 2, 3))
```

### 9. `mfield(dict_shape3d, point)`
Calculates the magnetic fields of wires at a point.

- **Parameters**:
  - `dict_shape3d`: A dictionary to describe the wire objects.
                    dict_shape3d = {
                                        { 'curve_name': {'curve': shap3d object,
                                                         'current': electric current through the curve
                                                        },
                                          ...
                                        }
                                        ...
                                    }
 - `point`: the destination point where you want to calculate the magnetic field, (x, y, z): tuple, list or 1darray.
- **Returns**:
  - The tuple of the magnetic fields of objects contributing to the point. The first element is the dictionary of the magnetic field of each wire. The second element is the total magnetic field of all the wires.

    ({'curve_name': magnetic field at point,
       ...
     },
     total magnetic field contributed by all wires
    )

```python
mfield(dict_shape3d, (1, 2, 3))
```

### 10. `mpotential(dict_shape3d, point)`
Calculates the magnetic vector potentials of wires at a point.

- **Parameters**:
  - `dict_shape3d`: A dictionary to describe the wire objects.
                    dict_shape3d = {
                                        { 'curve_name': {'curve': shap3d object,
                                                         'current': electric current through the curve
                                                        },
                                          ...
                                        }
                                        ...
                                    }
 - `point`: the destination point where you want to calculate the magnetic vector potential, (x, y, z): tuple, list or 1darray.
- **Returns**:
  - The tuple of the magnetic fields of objects contributing to the point. The first element is the dictionary of the magnetic field of each wire. The second element is the total magnetic field of all the wires.

    ({'curve_name': magnetic vector potential at point,
       ...
     },
     total magnetic vector potential contributed by all wires
    )

```python
mpotential(dict_shape3d, (1, 2, 3))
```

### 11. `mutual_inductance(shape3d_1, shape3d_2)`
Calculates the mutual inductance of two wires.

- **Parameters**:
  - `shape3d_1`: Wire 1 in Shape3d
  - `shape3d_2`: Wire 2 in Shape3d
- **Returns**:
  A floating number of the mutual inductance of the wires.

```python
mutual_inductance(shape3d_1, shape3d_2)
```

### 12. `self_inductance(shape3d)`
Calculates the self inductance of a wire.

- **Parameters**:
  - `shape3d`: Wire in Shape3d
- **Returns**:
  A floating number of the self inductance of the wire.

```python
self_inductance(shape3d)
```

For any concerns, Please contact: gang.wang@greenville.edu
