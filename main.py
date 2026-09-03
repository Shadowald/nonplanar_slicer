import numpy as np  # Used throughout for arrays
from stl import mesh  # Used only to import mesh data
import sys  # Used to increase recursive limit
import vg  # Used only to determine adjacency
from copy import deepcopy  # Used to create layers
from geomdl import BSpline  # Used to create layers and tool paths?
from geomdl import utilities  # Used to create layers and tool paths?
from geomdl import tessellate  # Used to tessellate the evaluated surfaces from set control points
import math  # Used to calculate angles between 2 lines

# TODO: Remove these imports when code is finished; used only for testing
import time  # Used for timing the code
from matplotlib import pyplot  # Used to plot the facets to confirm algorithms work
from mpl_toolkits import mplot3d  # Used to plot the facets to confirm algorithms work

EDGE_ANGLE = 80  # The value that determines if a facet is part of a side surface; units degrees
PERMITTED_LAYER_HEIGHTS = [.1, .2, .3, .4, .5]  # Available layer height values; units mm
NURBS_DEGREES = 4  # The degree value used to calculate the NURBS surface in layer generation
MAX_CTRL_PTS = 5  # The maximum distance between control points; units mm

# TODO: Remove this variable when code is finished; only to save time in testing
LAYER_HEIGHT = .1  # Layer height for slices; unit mm


# This class offers more flexibility in accessing facet data such as specific sides and adjacency information
class Facet:
    def __init__(self, facet, normal):
        self.facet = facet
        self.normal = normal
        self.edge = []
        self.corner = None
        self.adjacent = {'01': None, '02': None, '12': None}
        self.sides = {'01': facet[:2], '02': np.array([facet[0], facet[-1]]), '12': facet[1:3]}

    def set_adjacency(self, key, value):
        if key[0] > key[1]:
            temp = list(key)
            temp[0], temp[1] = temp[1], temp[0]
            key = ''.join(temp)
        self.adjacent[key] = value

    def get_adjacency(self, key=None):
        if key is None:
            return self.adjacent.values()
        else:
            return self.adjacent[key]

    def set_edge(self, value):
        self.edge.append(value)

    def set_corner(self, value):
        self.corner = value


# This class offers more flexibility in accessing line information
class Edge:
    def __init__(self, side, f_index):
        self.side = side
        self.f_index = f_index
        self.c_index = None
        self.corner = False
        self.adjacent = {0: None, 1: None}
        self.angles = {0: None, 1: None}
        self.relation = {0: None, 1: None}

    def set_adjacency(self, key, value):
        self.adjacent[key] = value

    def get_adjacency(self, key=None):
        if key is None:
            return self.adjacent.values()
        else:
            return self.adjacent[key]

    def set_angles(self, key, value):
        self.angles[key] = value

    def set_relations(self, key, value):
        self.relation[key] = value
        if value == 'corner':
            self.corner = True

    def get_relation(self):
        if self.relation[0] == 'straight' and self.relation[1] == 'straight':
            return 'straight'

        elif self.corner:
            return 'corner'

        else:
            return 'curved'


def main():
    sys.setrecursionlimit(1500)
    is_sorted_dict, facets, stl_mesh = initialize()
    slice_the_model(is_sorted_dict, facets, stl_mesh)

    """
    Simplest decomposition for a 5-axis slicer:

    generate_tool_paths()
    generate_gcode()
    Export G-code
    """


def graph_layer(stl_mesh, layer):
    """
    Graphs a list of facets into a 3D space. Initial idea was to just show a single layer with this function,
    but can be used with any single list of facets.

    :param stl_mesh: Mesh data from STL file
    :param layer: A list of facets to be graphed
    :return:
    """
    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)

    # Load the STL files and add the vectors to the plot
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(layer))

    # Auto scale to the mesh size
    scale = stl_mesh.points.flatten('F')
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def graph_edge(stl_mesh, edge):
    """
    Graphs a list of lines into a 3D space to show proper placement of toolpaths.

    :param stl_mesh: Mesh data from STL file
    :param edge: A list of lines to be graphed
    :return:
    """
    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)

    # Load the STL files and add the vectors to the plot
    axes.add_collection3d(mplot3d.art3d.Line3DCollection(edge))

    # Auto scale to the mesh size
    scale = stl_mesh.points.flatten('F')
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def graph_layers(stl_mesh, layers):
    """
    Graphs a list of a list of facets into a 3D space. This similar to the graph_layer function,
    except that it is designed to graph multiple layers as opposed to one. This requires the list
    to be passed in to be populated with lists of facets, like how the graph_layer function was
    intended to be used.

    :param stl_mesh: Mesh data from STL file
    :param layers: A list of facets to be graphed
    :return:
    """
    # Create a new plot
    figure = pyplot.figure()
    # axes = figure.add_subplot(1, 1, 1, projection='3d')
    axes = mplot3d.Axes3D(figure)

    color = ['b', 'r', 'g', 'm', 'c', 'y']

    # Load the STL files and add the vectors to the plot
    for clr, layer in enumerate(layers):
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(layer))
        # axes.add_collection3d(mplot3d.art3d.Poly3DCollection(layer, facecolor=color[clr]))

    """axes = figure.add_subplot(1, 1, 2, projection='3d')
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(stl_mesh.vectors))"""

    # Auto scale to the mesh size
    scale = stl_mesh.points.flatten('F')
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def read_stl():
    """
    This function reads in the data from an stl file of the same name inputted by the user.
    It is important that the file to be read in is already within the SlicerMain folder
    before running the program.

    :return: the mesh of the uploaded stl file
    """

    # TODO: Reinstate lines below once program is finished
    """
    while True:
        # Tries to open a file, based on a value given by the user
        try:
            file_name = input("Please enter the name of the stl file you wish to open with the file extension:  ")
            stl_mesh = mesh.Mesh.from_file(file_name)
        # If no such file is found, exception is thrown and loop restarts
        except FileNotFoundError:
            print("Sorry, no such file was found.")
            print("Please check that you entered the name correctly and/or added .STL at the end")
            # better try again... Return to the start of the loop
            continue
        else:
            # Got the correct format
            # we're ready to exit the loop.
            break
    """

    stl_mesh = mesh.Mesh.from_file('Plane1.STL')
    return stl_mesh


def get_layer_height():
    """
    Asks the user for their desired layer height for their free form surface to be printed with.
    If value given is not of the type float or in the permitted value list, an exception is
    thrown and user is asked again.

    :return: The layer height value given by the user
    """

    while True:
        # Tries to get a float value from user
        try:
            layer_height = float(input("Please enter your desired layer height in mm: "))
        # If value is not a float, then exception is thrown and loop restarts
        except ValueError:
            print("Sorry, I didn't understand that.")
            continue

        # If value is a float, but one of the accepted values, then user is informed of allowed values and loop restarts
        if layer_height not in PERMITTED_LAYER_HEIGHTS:
            print("That is not an accepted layer height. Accepted heights are:")
            print(*PERMITTED_LAYER_HEIGHTS, sep=', ')
            continue
        else:
            # layer height was successfully parsed, and we're happy with its value.
            # we're ready to exit the loop.
            break

    return layer_height


def initialize():
    """
    Initializes the various objects used throughout the program, while also saving some of the memory space.

    :return:
    """

    stl_mesh = read_stl()
    # TODO: Reinstate later, only to save time in testing
    # layer_height = get_layer_height()

    is_sorted_dict = make_dict(stl_mesh)
    facets = []  # A list that holds all the objects of type Facet

    # Populates list facets with all the facets and their relevant data
    for i, f1 in enumerate(stl_mesh.vectors):
        facets.append(Facet(f1, stl_mesh.normals[i]))

    # TODO: Add layer_height to the return list and remove stl_mesh
    return is_sorted_dict, facets, stl_mesh


def get_curves(is_sorted_dict):
    curved_surfaces = []
    print("Separating curved surfaces from side surfaces.")
    sort_via_dict(is_sorted_dict, curved_surfaces)
    print("Separation complete.")
    return curved_surfaces


def slice_the_model(is_sorted_dict, facets, stl_mesh):
    """
    Slices the STL mesh into curved layers.

    :param is_sorted_dict: A dictionary that keeps track of whether the facets have been sorted
    :param facets: A list with facet metadata
    :return: Layer slices
    """

    side_surfaces = detect_side_surfaces(is_sorted_dict, facets)
    curved_surfaces = get_curves(is_sorted_dict)
    set_edges(facets, curved_surfaces, side_surfaces)
    top_layer, bottom_layer = separate_layers(curved_surfaces, facets, is_sorted_dict)

    layer_slices = create_layers(top_layer, bottom_layer, facets, side_surfaces, stl_mesh)

    return layer_slices


def detect_side_surfaces(is_sorted_dict, facets):
    """
    From the STL mesh, determines which facet vectors make up the
    4 side surfaces and which make up the 2 curved surfaces.

    :param is_sorted_dict: A dictionary that keeps track of whether the facets have been sorted
    :param facets: A list with facet metadata
    :return: the list of objects type Facet
    """

    side_surfaces = []  # List that will hold all side surface facets

    # Another list that was used to extract data and improve global variable values
    normals_angle = []

    # Print statement to let the user know program is running
    print("Detecting side surfaces. Please wait.")
    for i in range(0, len(facets)):
        f1 = facets[i].facet  # Holds the facet data at the given index

        # Nested for loop is intentionally made to start one after the outer for loop
        # to save time by not checking earlier facets
        for j in range(i + 1, len(facets)):
            f2 = facets[j].facet  # Holds the facet data at the given index

            # Prevents duplicate bug
            if (f1 == f2).all():
                continue

            # Skips checking adjacency if all sides have had adjacency set
            if all(facets[i].get_adjacency()):
                break

            is_adjacent, side_str1, side_str2 = check_adjacency(f1, f2)
            if is_adjacent:
                facets[i].set_adjacency(side_str1, j)
                facets[j].set_adjacency(side_str2, i)

                # Determines angle between the 2 adjacent normals in units degrees
                angle = vg.angle(facets[i].normal, facets[j].normal, units='deg')

                normals_angle.append(angle)

                # if angle is greater than EDGE_ANGLE, than outer loop facet is part of the side surfaces
                # Second conditional prevents duplication bug with corner facets
                if angle > EDGE_ANGLE and not is_sorted_dict[i]:
                    side_surfaces.append(i)
                    is_sorted_dict[i] = True

    print("Side surfaces found.")
    """
    Writes the list of angle differences between adjacent normals to a txt file.
    This data was used to help determine appropriate values for ZERO and MAX_ANGLE-DIFFERENCE.

    with open('angle list.txt', 'w') as filehandle:
        for listitem in normals_angle:
            filehandle.write('%s\n' % listitem)
    """

    return side_surfaces


def check_adjacency(f1, f2):
    """
    Determines if two facets are adjacent to each other. Was coded to work with any tesselation shape.

    :param f1: first facet
    :param f2: second facet
    :return:is_adjacent: a bool value, True if adjacent, False if not
            side_str1, side_str2: a string that represents which sides
            of the facet are adjacent based on index values for facet pts
            ie, '01', '02', '12'
    """

    count_matching_pts = 0  # Counts how many vertices have matched
    is_adjacent = False  # A bool representing whether facets are adjacent
    side_str1 = ""  # A string that represents which side of f1 is adjacent to f2
    side_str2 = ""  # A string that represents which side of f2 is adjacent to f1

    # Nested for loop which compares each point of one facet with each point of another facet, looking for matching pts
    for i in range(len(f1)):
        for j in range(len(f2)):

            # If all values of a point match with another, the counter is incremented by one
            # and index values are added to corresponding strings
            if (f1[i] == f2[j]).all():
                count_matching_pts += 1
                side_str1 += str(i)
                side_str2 += str(j)

            # Once 2 matching points are found, the facets are concluded to be adjacent
            if count_matching_pts == 2:
                is_adjacent = True

                return is_adjacent, side_str1, side_str2

    return is_adjacent, side_str1, side_str2


def set_edges(facets, curved_surfaces, side_surfaces):
    # Gets the index in layer
    for f1 in curved_surfaces:

        # Using the index, get the key to the adjacent dictionary associated with that facet
        for key in facets[f1].adjacent:

            # Gets the index in side_surfaces
            for f2 in side_surfaces:

                # Check if the adjacent dictionary value is the same as the index value from side_surfaces
                if facets[f1].adjacent[key] == f2:
                    # Sets side data as an edge
                    facets[f1].set_edge(facets[f1].sides[key])
                    break


def separate_layers(curved_surfaces, facets, is_sorted_dict):
    """
    Separates the top and bottom surfaces into their own lists

    :param curved_surfaces: The list of facets that makes up the curved surfaces
    :param facets: List containing relevant facet data
    :param is_sorted_dict: The dictionary of index keys and boolean values of whether or not it was sorted
    :return: Two lists, one containing all the facets for the top surface, and one for the bottom surface
    """

    bottom_lyr = [curved_surfaces[-1]]  # A list to hold all facets that make up the top most surface of the STL

    # Changes the value in the is_sorted_dict to True
    is_sorted_dict[curved_surfaces[-1]] = True

    print("Separating top layer from bottom layer.")
    sort_via_adjacency(curved_surfaces[-1], facets, is_sorted_dict, bottom_lyr)

    curved_surfaces = [x for x in curved_surfaces if x not in bottom_lyr]

    top_lyr = [curved_surfaces[-1]]

    # Changes the value in the is_sorted_dict to True
    is_sorted_dict[curved_surfaces[-1]] = True

    sort_via_adjacency(curved_surfaces[-1], facets, is_sorted_dict, top_lyr)

    curved_surfaces = [x for x in curved_surfaces if x not in top_lyr]

    # Using the logic that only the bottom surface facets remain 'unsorted',
    # they are now officially sorted into bottom_layer.
    # bottom_lyr = []  # A list to hold all facets that make up the bottom most surface of the STL
    # sort_via_dict(is_sorted_dict, bottom_lyr, say_sorted=True)
    print("Layers have been separated.")

    return top_lyr, bottom_lyr


def sort_via_adjacency(top, facets, is_sorted_dict, top_lyr):
    """
    A recursive function that takes in a facet index value, among other needed lists and dicts, that checks for
    the adjacency values stored. If the value has not been listed as sorted, then it is appended to top_lyr,
    its corresponding value in is_sorted_dict is changed to True, and the function is called passing in newly
    sorted value as the new facet index value.

    :param top: index of facet to get metadata
    :param facets: List containing relevant facet data
    :param is_sorted_dict: The dictionary of index keys and boolean values of whether or not it was sorted
    :param top_lyr: list containing indexes of facets that are apart of the top surface
    :return: No return value due to passing in by reference
    """

    # Gets the key to the adjacent dictionary associated with facet at index top
    for key in facets[top].adjacent:

        # Stores dictionary value using key in variable index, which is an index associated with a facet
        index = facets[top].get_adjacency(key)

        # If facet associated with index has yet to be sorted, add index to top_lyr,
        # change is_sorted_dict value to True, and recursively call this function
        # again with index as the first given parameter.
        if not is_sorted_dict[index]:
            top_lyr.append(index)
            is_sorted_dict[index] = True
            sort_via_adjacency(index, facets, is_sorted_dict, top_lyr)


def make_dict(stl_mesh):
    """
    Makes a dictionary in which the keys are the indices of the facets in the
    stl mesh and the values are a bool set to False.

    :param stl_mesh: Mesh data from STL file
    :return: The dictionary of indices
    """
    dict_of_ids = {}
    for i in range(len(stl_mesh.vectors)):
        dict_of_ids[i] = False
    return dict_of_ids


def get_top_facet(curved_surface, facets):
    """
    Finds the top most facet, in terms of z value, in the stl mesh.

    :param facets: List containing relevant facet data
    :param curved_surface: A list of sll facets that ar+e apart of the curved surfaces
    :return: The top most facet and its index value in the stl mesh
    """
    top_facet = np.zeros((3, 3))  # A variable to hold the highest facet in terms of z-value

    # Checks only the z values within the facet. When a facet is found
    # with at least one z-point higher than the current top_facet,
    # that facet becomes the new top_facet
    for item in curved_surface:
        temp_facet = facets[item].facet
        if (temp_facet[0:3, 2] > top_facet[0:3, 2]).any():
            top_facet = facets[item].facet
            index = item
    return top_facet


def sort_via_dict(is_sorted_dict, list_to_sort, say_sorted=False):
    """
    Sorts facets based on what has already been sorted according to the is_sorted dictionary.

    :param is_sorted_dict: The dictionary of index keys and boolean values of whether or not it was sorted
    :param list_to_sort: The list in which the facets will be sorted into
    :param say_sorted: A boolean value which says whether to change the dictionary to True during sorting.
           Default value is False
    :return: No return necessary as dictionaries and lists are mutable objects
    """

    # Iterates through all keys in the dictionary
    for key in is_sorted_dict:

        # If the value attached to that key is False,
        # the facet associated with that key is added to the list passed through

        if not is_sorted_dict[key]:
            list_to_sort.append(key)

            # If say_sorted is passed in as True, changes dictionary value to True for current key
            if say_sorted:
                is_sorted_dict[key] = True


def create_layers(top, bottom, facets, side_surfaces, stl_mesh):
    """
    Creates the layers to generate the tool paths onto

    :param bottom: list of indexes of all facets comprising of the bottom layer
    :param top: list of indexes of all facets comprising of the top layer
    :param facets: List containing relevant facet data
    :param side_surfaces: list of indexes of all facets comprising of the side surfaces
    :return: A list of the layers of the model
    """
    layer_slices = []  # List to hold all layers of facets
    # layer = []
    temp = []

    print("Creating layer slices.")
    # Calculates the total number of layers to be generated
    range_num = calc_num_layers(top, bottom, side_surfaces, facets)
    bttm = []
    for f in bottom:
        bttm.append(Facet(facets[f].facet, facets[f].normal))

        temp.append(facets[f].facet)

    """# Gets the edges of layer, then sets the control points
    edges = find_edges(facets, bottom, side_surfaces)
    control_pts, num_x, num_y = get_control_pts(facets, bottom, edges)

    # Evaluates the layer using the control points
    surf = evaluate_surface(control_pts)

    for face in surf.faces:
        vector = []
        for vertex in face.vertices:
            vert = np.array(vertex)
            vector.append(vert)
        layer.append(vector)

    graph_layer(stl_mesh, layer)"""

    # TODO: Fix code so that the new layers made work with the code implemented for already existing surfaces
    # Layers are created by going over the total number of layers minus 2 (for the top and bottom)
    edges = find_edges_id(facets, bottom, side_surfaces)
    control_pts, control_pts_id, num_x, num_y = get_control_pts(facets, bottom, edges)
    layer = []

    for i in range(range_num - 2):
        # First time through uses bottom layer, every other time uses the previous made layer
        """if i == 0:
            # layer = create_layer(facets, bottom)
            graph_layer(stl_mesh, temp)
            layer = create_layer(bttm)
            temp.clear()
            for f in layer:
                temp.append(f.facet)
            graph_layer(stl_mesh, temp)
            # edges = find_edges_math(facets, layer, side_surfaces)
            # graph_edge(stl_mesh, edges)
        else:
            # layer = create_layer(facets, layer_slices[i-1])
            layer = create_layer(layer_slices[i - 1])

            # Gets the edges of layer, then sets the control points
            # edges = find_edges_math(facets, layer, side_surfaces)
            # graph_edge(stl_mesh, edges)"""

        # control_pts, num_x, num_y = get_control_pts(facets, layer, edges)
        control_pts = translate_ctrl_pts(control_pts, control_pts_id, facets)

        # Evaluates the layer using the control points
        surf = evaluate_surface(control_pts)
        graph_layer(stl_mesh, surf)

        # Layer is cleared and then repopulated using facets from the surface
        layer.clear()
        for face in surf.faces:
            vector = []
            for vertex in face.vertices:
                vert = np.array(vertex)
                vector.append(vert)
            layer.append(vector)

        # Finally layer is deepcopied into layer_slices, to avoid list of duplicates bug
        layer_slices.append(deepcopy(layer))

    # TODO: make it so the code below adds a list of facets and not indexes to the list
    # layer_slices.insert(0, bottom)
    # layer_slices.append(top)

    # TODO: Determine if code below should be integrated or not
    """
    Makes the final layer equal to the difference between the last layer generated 
    and the top layer from stl, if the the difference is greater than or equal to
    half of the layer height.    

    index = get_top_facet(bottom, stl_mesh)
    top_bottom_facet = temp_layer[bottom.index(index)]
    diff = min(top_facet[:, 2] - top_bottom_facet[:, 2])
    if diff >= (LAYER_HEIGHT / 2):
        temp_layer = []
        for facet in bottom:
            temp_facet = deepcopy(stl_mesh.vectors[facet])
            for point in temp_facet:
                point[2] = point[2] + (LAYER_HEIGHT * i) + diff
            temp_layer.append(temp_facet)
        layer_slices.append(temp_layer)
    """

    return layer_slices


def translate_ctrl_pts(ctrl_pts, ctrl_pts_ids, facets):
    new_ctrl_pts = []

    for i in range(len(ctrl_pts)):
        temp = []
        for j in range(len(ctrl_pts[i])):
            increment = facets[ctrl_pts_ids[i, j]].normal * LAYER_HEIGHT
            new_pt = deepcopy(ctrl_pts[i, j])
            new_pt += increment
            temp.append(new_pt)
        new_ctrl_pts.append(deepcopy(temp))

    return new_ctrl_pts


def calc_num_layers(top, bottom, side_surfaces, facets):
    """
    Calculates how many layers will need to be generated. This is done by using the pythagoras theorem
    to find the distance between the corners closest to the origin and divides this distance by the set layer height.

    :param top: A list of indexes representing all facets on the top surface
    :param bottom: A list of indexes representing all facets on the bottom surface
    :param side_surfaces: A list of indexes representing all facets on the side surfaces
    :param facets: List containing relevant facet data
    :return: A int value that represent the total number of layers needed
    """
    # Gets the edges of the top and bottom layers
    e1 = find_edges_id(facets, top, side_surfaces)
    e2 = find_edges_id(facets, bottom, side_surfaces)

    # Uses the edges to find the corners in their respective layers
    c1, d1 = find_corners(e1)
    c2, d2 = find_corners(e2)

    # Stores the first corner in each list for readability.
    # Any corner could be used so long its the same position in both lists.
    p1 = c1[0]
    p2 = c2[0]

    # Calculates the difference between each coordinate and squares the difference
    x = math.pow(p2[0] - p1[0], 2)
    y = math.pow(p2[1] - p1[1], 2)
    z = math.pow(p2[2] - p1[2], 2)

    # Calculates the distance between the points using the pythagoras theorem and rounds the answer
    diff = round(math.sqrt(x + y + z))

    # Total number of layers calculated by dividing diff by LAYER_HEIGHT
    return int(diff // LAYER_HEIGHT)


# def create_layer(facets, layer):
def create_layer(layer):
    """
    Creates a new layer by incrementing by the LAYER_HEIGHT across the normal vectors of every facet

    :param facets: List of all facets
    :param layer: List of of all facets in layer
    :return: The new layer
    """

    new_layer = []  # List to hold all the new facets

    # TODO: Modify code to work for all layers. Currently can make the first layer after the bottom,
    #       but will not work for all subsequent layers.
    # TODO: Modify code so increment is equal to average of facet normals*LAYER_HEIGHT associated with the point.
    # Iterates over each facet and uses their normal vector multiplied by the LAYER_HEIGHT to create increment vector
    for f in layer:
        # increment = facets[f].normal * LAYER_HEIGHT
        increment = f.normal * LAYER_HEIGHT

        # Deepcopies the facet to prevent retroactive alteration bug
        # new_facet = deepcopy(facets[f].facet)
        new_facet = deepcopy(f)

        # Every vertex has the increment vector added to it
        for vert in new_facet.facet:
            vert += increment

        # Adds the new facet to the layer
        new_layer.append(new_facet)

    return new_layer


def find_corners(edges):
    """
    Finds the corners of the of the layer using the edge data

    :param edges: A list of line segments that make up the edge of the layer.
    :return: A list of the corner points.
    """
    corners = []  # A list that will hold the corner points
    c_index = []  # Holds the index values of the corners in pairs
    angles = []

    for i in range(len(edges)):
        e1 = edges[i].side
        count = 0  # Counts the number of matched points to speed up run time

        for j in range(i + 1, len(edges)):
            e2 = edges[j].side
            # A boolean of whether the 2 segments are adjacent and the shared point
            is_adjacent, point, in1, in2 = line_adjacency(e1, e2)

            if is_adjacent:
                count += 1

                edges[i].set_adjacency(in1, j)
                edges[j].set_adjacency(in2, i)

                # Calculate the vector representations of the line segments, used to calculate angle
                v1, v2 = line_2_vector(e1, e2)

                # Calculate the angle between the line segment, range 0 - 180 degrees
                angle = calc_angle(v1, v2)
                angles.append(angle)

                edges[i].set_angles(in1, angle)
                edges[j].set_angles(in2, angle)

                # Tolerance of 10 degrees was used to differentiate corners from curve approximations
                if 10 < angle < 170:
                    corners.append(deepcopy(point))
                    # corners.append(point)

                    # Appends indexes after sorting them with the greater difference in the x-direction first
                    c_index.append(greater_x_diff(e1, e2, i, j))

                    edges[i].set_relations(in1, 'corner')
                    edges[j].set_relations(in2, 'corner')

                elif angle == 0 or angle == 180:
                    edges[i].set_relations(in1, 'straight')
                    edges[j].set_relations(in2, 'straight')

                else:
                    edges[i].set_relations(in1, 'curved')
                    edges[j].set_relations(in2, 'curved')

            # Speeds up run time
            if count == 2:
                break

    # Sort the corners so that the bottom left hand corner from the XY perspective is in the 0th position,
    # with the subsequent positions placed clockwise around the surface.
    c_id_dict = sort_corners(corners, c_index)

    """with open('angle list2.txt', 'w') as filehandle:
        for listitem in angles:
            filehandle.write('%s\n' % listitem)"""

    return corners, c_id_dict


def line_adjacency(e1, e2):
    """
        Determines if two edge lines are adjacent to each other.

        :param e1: first edge line
        :param e2: second edge line
        :return: is_adjacent: a bool value, True if adjacent, False if not
                 The adjacent point, if it exists, is the second return parameter
        """

    is_adjacent = False  # A bool representing whether facets are adjacent

    # Nested for loop which compares each point of one facet with each point of another facet, looking for matching pts
    for i in range(len(e1)):
        for j in range(len(e2)):

            # If all values of a point match with another, the counter is incremented by one
            # and index values are added to corresponding strings
            if (e1[i] == e2[j]).all():
                is_adjacent = True
                return is_adjacent, e1[i], i, j

    # In the event that the lines are not adjacent, a value of None is returned
    return is_adjacent, None, None, None


def line_2_vector(e1, e2):
    """
    Calculates the vector representations of the given line segments via linear algebra.

    :param e1: first line segment or edge
    :param e2: second line segment or edge
    :return: v1: the vector representation of the first edge
             v2: the vector representation of the second edge
    """

    # Vector calculated by subtracting the first coordinate from the respective second coordinate
    v1 = [e1[1][0] - e1[0][0], e1[1][1] - e1[0][1], e1[1][2] - e1[0][2]]
    v2 = [e2[1][0] - e2[0][0], e2[1][1] - e2[0][1], e2[1][2] - e2[0][2]]

    return v1, v2


def calc_angle(v1, v2):
    """
    Calculate the angle between 2 line segments via their vector representations with the use of the dot product and
    arc cosine.

    :param v1: First vector
    :param v2: Second vector
    :return: The angle between the lines in units degrees
    """

    # Separates the vector into its x-, y-, and z-coordinates
    x1, y1, z1 = v1
    x2, y2, z2 = v2

    # The dot product calculation
    inner_product = x1 * x2 + y1 * y2 + z1 * z2

    # Calculate the length of the vector
    len1 = math.hypot(x1, y1, z1)
    len2 = math.hypot(x2, y2, z2)

    # Calculate the angle in radians
    rads = math.acos(round(inner_product / (len1 * len2), 3))

    # Return the angle in degrees
    return abs(math.degrees(rads))


def greater_x_diff(e1, e2, i, j):
    """
    Determines which edge has the greatest difference in the x-coordinate.

    :param e1: First edge
    :param e2: Second edge
    :param i: First index
    :param j: Second index
    :return: List of the indexes ordered
    """

    diff1 = abs(e1[0][0] - e1[1][0])
    diff2 = abs(e2[0][0] - e2[1][0])

    if diff1 < diff2:
        return [j, i]
    else:
        return [i, j]


def sort_corners(corners, c_index):
    """
    Sorts the corner points so that the bottom, left-hand corner is in the 0th position of the list,
    with the 1st position holding the top, left-hand corner, and continuing in a clockwise fashion.

    :param c_index: A list that holds the index values of the corners
    :param corners: The list unsort corner points
    :return: no return value due to passing in by reference
    """

    temp = np.array(corners)  # An array to hold the corner point data, used for slicing purposes
    c_id_dict = {}  # A dictionary to hold the indexes of the edges that belong to each corner

    # Sorts the x- and y-coordinates from smallest to largest and  stores their index value in a list,
    # with the index of smallest coordinate in position 0
    srtd_x_coords = smallest_2_largest(temp[:, 0:1])
    srtd_y_coords = smallest_2_largest(temp[:, 1:2])

    # Using the logic that the smallest values are in position 0 and 1, and the largest in 3 and 4,
    # sorts the points in the corners list
    for i in range(len(srtd_x_coords)):
        for j in range(len(srtd_y_coords)):
            if srtd_x_coords[i] == srtd_y_coords[j]:
                index = srtd_x_coords[i]
                if i < 2 and j < 2:
                    corners[0] = temp[index]
                    c_id_dict[0] = c_index[index]
                elif i < 2 and j > 1:
                    corners[1] = temp[index]
                    c_id_dict[1] = c_index[index]
                elif i > 1 and j > 1:
                    corners[2] = temp[index]
                    c_id_dict[2] = c_index[index]
                elif i > 1 and j < 2:
                    corners[3] = temp[index]
                    c_id_dict[3] = c_index[index]

    return c_id_dict


def smallest_2_largest(coordinates):
    """
    Sorts an array of 4 coordinates from smallest to largest.

    :param coordinates: A array of 4 coordinates
    :return: A list of the indexes sorted from smallest to largest based on the values they represent within the array
    """

    largest = None  # Holds the largest value
    lowest = None  # Holds the lowest value
    largest2 = None  # Holds the 2nd largest value
    lowest2 = None  # Holds the 2nd lowest value

    # Sorts the coordinate values by assigning them to the variables made above
    for i in range(len(coordinates)):
        if largest is None or coordinates[i] > largest:
            largest2 = largest
            largest = coordinates[i]
        elif largest2 is None or largest2 < coordinates[i]:
            largest2 = coordinates[i]
        if lowest is None or coordinates[i] < lowest:
            lowest2 = lowest
            lowest = coordinates[i]
        elif lowest2 is None or lowest2 > coordinates[i]:
            lowest2 = coordinates[i]

    # Changes the array to a list, done to use .index() function
    corr_list = coordinates.tolist()

    # Assigns the index of the found values to variables. If-else statements are used to avoid duplication bug
    i1 = corr_list.index(lowest)

    if lowest == lowest2:
        for i in range(len(corr_list)):
            if corr_list[i] == lowest2:
                i2 = i
    else:
        i2 = corr_list.index(lowest2)

    i3 = corr_list.index(largest2)

    if largest == largest2:
        for i in range(len(corr_list)):
            if corr_list[i] == largest:
                i4 = i
    else:
        i4 = corr_list.index(largest)

    return [i1, i2, i3, i4]


def set_edge_ctrl_pts(edges, num_cpts, coord, pos, end):
    """
    Sets the the control points along the edge of the layer.

    :param edges: List of type class Edge
    :param num_cpts: Calculated number of control points in given direction
    :param coord: Integer to represent main coordinate, 0 for x and 1 for y
    :param pos: Position of current control point being set
    :param end: End point of where to set control points
    :return: List of the layer edge control points
    """

    edge_pairs = []  # List that holds a pairs of points
    edge_pairs_id = []

    # Spacing between control points is calculated
    spacing = (end - pos) / num_cpts

    # while loops runs the length of the layer in the  coord direction
    while pos < end:
        edge_pair = []  # Holds the pair of points to be set
        edge_pair_id = []

        # x is already incremented as corners are already stored
        pos += spacing

        # Used to prevent the distance between final point and the corner from being unnecessarily small
        if (pos + (spacing * 0.2)) > end:
            break

        for edge in edges:

            # Converts the list of arrays to a list of lists in order to properly use sort()
            line = edge.side.tolist()
            line.sort(key=lambda a: a[coord])

            # Converted to an array of arrays for computational purposes
            line = np.array(line)

            # When the x-coord lies between 2 points the equation for a line in 3D space is used to
            # calculate the appropriate y- and z-coords, after calculating the needed variables.
            # The coords are then combined into a list and appended to edge_pair
            # Equ. of a line in 3D space: [x, y, z] = p0 + t * r
            if line[0, coord] <= pos <= line[1, coord]:
                r = line[1] - line[0]
                t = (pos - line[0, coord]) / r[coord]
                z = line[0, 2] + (r[2] * t)

                # To calculate the correct coordinate depending on whether direction is in x (0) or y (1)
                if coord == 0:
                    y = line[0, 1] + (r[1] * t)
                    edge_pair.append([pos, y, z])
                    edge_pair_id.append(edge.f_index)
                elif coord == 1:
                    x = line[0, 0] + (r[0] * t)
                    edge_pair.append([x, pos, z])
                    edge_pair_id.append(edge.f_index)

            # Speeds up run time
            if len(edge_pair) == 2:
                pass  # break

        # Sorts the pair so that the lesser coord is first
        if coord == 0:
            temp = edge_pair[0]
            edge_pair.sort(key=lambda a: a[1])
            if temp != edge_pair[0]:
                edge_pair_id.reverse()
        if coord == 1:
            temp = edge_pair[0]
            edge_pair.sort(key=lambda a: a[0])
            if temp != edge_pair[0]:
                edge_pair_id.reverse()
        # edge_pair = np.array(edge_pair)

        # Adds a deepcopy of the edge pair to the dict with the current index value, which is then increased by 1
        edge_pairs.append(deepcopy(edge_pair))
        edge_pairs_id.append(deepcopy(edge_pair_id))

    return edge_pairs, edge_pairs_id


def is_side_straight(edges, index, pos, end, coord):
    """
    Determines if the side of a quadrilateral is curved or straight through a recursive method.

    :param edges: List of class type Edge
    :param index: The index of the Edge in edges being checked.
    :param pos: The current position on the side.
    :param end: The end point, AKA the next corner.
    :param coord: Integer to represent main coordinate, 0 for x and 1 for y
    :return: A list of True and False values that indicate if the side is straight
    """

    straight = True  # straight is initially set as True

    # While loop ends if straight is changed to False or pos is equal to end
    while straight and (pos != end).any():

        # Straightness is determined by the main coord of pos being nearly the same as the end coord
        # If straight, updates necessary values and calls itself, else updates straight value and exists recursive loop
        if abs(end - pos)[coord] < 1:

            line = edges[index]  # Holds line data from the edges to make code cleaner

            # Goes through the points in line to find the one that is not pos.
            # Then uses this to update the the index and pos value, then recursively calls itself.
            # Break statement is used to prevent infinite recursive loop.
            for i in range(len(line.side)):
                if (line.side[i] != pos).all():
                    index = line.adjacent[i]
                    pos = line.side[i]
                    straight, pos = is_side_straight(edges, index, pos, end, coord)
                    break

        else:
            straight = False

    return straight, pos


def get_side_relations(edges, c_id_dict, corners):
    """
    Determines if the sides of a quadrilateral are straight or curved.

    :param edges: List of class type Edge
    :param c_id_dict: Dict of indexes for corners.
    :param corners: A sorted list of the corners.
    :return: A list of boolean values that represent whether a side is straight
    """

    side_relations = []  # Holds the boolean values

    # Adds the first entry to the end to make it loop back for later code.
    corners.append(corners[0])
    c_id_dict[4] = c_id_dict[0]

    # Range has been hard-coded for quadrilaterals
    # TODO: Make code usable for all shapes
    for i in range(4):

        # Uses remainder to determine which side is being checked
        if i % 2 == 0:
            straight, ignore = is_side_straight(edges, c_id_dict[i][1], corners[i], corners[i + 1], 0)
        else:
            straight, ignore = is_side_straight(edges, c_id_dict[i][0], corners[i], corners[i + 1], 1)
        side_relations.append(straight)

    # Removes the position that was added on at beginning of the function
    corners.pop(4)
    c_id_dict.pop(4)

    return side_relations


def get_control_pts(facets, layer, edges):
    """
    Sets the control points on the layer and organizes them into a 2D grid pattern,
    ignoring that the points technically make it a 3D list

    :param facets: List of Class type Facet
    :param layer: List indexes of facets in layer
    :param edges: List of Class type Edges relating to layer
    :return: 3D list of control points. Control points are in a 2D grid pattern
    """

    # Find the corners of the layer
    corners, c_id_dict = find_corners(edges)

    # Using the corners, calculates the number of control points in both the x- and y-directions
    num_cpts_x, num_cpts_y = calc_num_pts(corners)

    # Determines if layer sides are straight or curved
    side_relations = get_side_relations(edges, c_id_dict, corners)

    # Setting the edge control points in pairs in both the x- and y-direction
    x_edge_pairs, x_edge_pairs_id = set_edge_ctrl_pts(edges, num_cpts_x, 0, corners[0][0], corners[3][0])
    y_edge_pairs, y_edge_pairs_id = set_edge_ctrl_pts(edges, num_cpts_y, 1, corners[0][1], corners[1][1])

    # Depending on the results of the side relations, determines in which direction to set the inner control points.
    # They are set along the curved sides, so that they run parallel to the straight sides.
    if not side_relations[1] or not side_relations[3]:
        inner_ctrl_pts, inner_ctrl_pts_id = set_inner_ctrl_pts(facets, layer, x_edge_pairs, num_cpts_y, 1, 0)
        inner_dirc = 1
    elif not side_relations[0] or not side_relations[2]:
        inner_ctrl_pts, inner_ctrl_pts_id = set_inner_ctrl_pts(facets, layer, y_edge_pairs, num_cpts_x, 0, 1)
        inner_dirc = 0
    elif len(y_edge_pairs) == 0 and len(x_edge_pairs) == 0:
        inner_ctrl_pts, inner_ctrl_pts_id = set_inner_ctrl_pts(facets, layer, y_edge_pairs, num_cpts_x, 0, 1)
        inner_dirc = 0
        # TODO: rewrite to work with NURBS library
        """inner_ctrl_pts = inner_pts_exception_case(facets, layer, spacing_x, spacing_y, num_cpts_x, num_cpts_y, corners)
        inner_dirc = None"""

    # Converts corner entries to list so that surface evaluation works correctly
    for item in corners:
        item.tolist()

    # Organizes the control points into a 2D grid pattern so that surface evaluation works
    control_pts, control_pts_id = organize_ctrl_pts(corners, c_id_dict, x_edge_pairs, x_edge_pairs_id, y_edge_pairs, y_edge_pairs_id, inner_ctrl_pts, inner_dirc, inner_ctrl_pts_id)

    """with open('ctrl pt index.txt', 'w') as filehandle:
        for listitem in pt_index:
            filehandle.write('%s\n' % listitem)"""

    return control_pts, control_pts_id, num_cpts_x, num_cpts_y


def organize_ctrl_pts(corners, c_ids, x_pairs, x_ids, y_pairs, y_ids,  inner_pts, dirc, inner_ids):
    """
    Organizes all the control points into a 2D grid in order for the geomdl library to evaluate surfaces.

    :param corners: List of corner points
    :param x_pairs: List of edge pair control points in the x direction
    :param y_pairs: List of edge pair control points in the y direction
    :param inner_pts: List of lists of the inner control points
    :param dirc: Variable to keep track of the direction in which the inner control points were set
    :return: A 2D grid layout of all the control points
    """

    ctrl_pts = []  # A list to hold all of the control points
    ctrl_pts_ids = []

    # Depending of the direction in which the inner points are set, determines how the grid is constructed.
    # It is constructed so that the points are added in the same direction as the inner points were set.
    if dirc == 0:
        outer_pairs = x_pairs
        outer_ids = x_ids
        end1 = 3
        inner_pairs = y_pairs
        inner_ids = y_ids
        start2 = 1

    elif dirc == 1:
        outer_pairs = y_pairs
        outer_ids = y_ids
        end1 = 1
        inner_pairs = x_pairs
        inner_ids = x_ids
        start2 = 3

    # Temp list is made to hold all points in the construction direction before being append to ctrl_pts
    # Construction always begins with the first point in the corner list and always end with the third point.
    temp = [corners[0].tolist()]
    temp_id = [c_ids[0]]
    for value in outer_pairs:
        temp.append(value[0])
    for id in outer_ids:
        temp_id.append(id[0])
    temp.append(corners[end1].tolist())
    temp_id.append(c_ids[end1])
    ctrl_pts.append(deepcopy(temp))
    ctrl_pts_ids.append(deepcopy(temp_id))

    # Adds the first point in the ith position of inner_pairs,
    # then all points in the ith position of inner_pts, and
    # finally the last point in the ith position of inner_pairs.
    # Reverse indexing, -1, is used to avoid error of pairs having more than 2 values.
    for i in range(len(inner_pairs)):
        temp.clear()
        temp_id.clear()
        temp.append(inner_pairs[i][0])
        temp_id.append(inner_ids[i][0])
        for pts in inner_pts[i]:
            temp.append(pts)
        for id in inner_ids[i]:
            temp_id.append(id)
        temp.append(inner_pairs[i][-1])
        temp_id.append(inner_ids[i][-1])
        ctrl_pts.append(deepcopy(temp))
        ctrl_pts_ids.append(deepcopy(temp_id))

    temp.clear()
    temp = [corners[start2].tolist()]
    temp_id = [c_ids[start2]]
    for value in outer_pairs:
        temp.append(value[-1])
    for id in outer_ids:
        temp_id.append(id[-1])
    temp.append(corners[2].tolist())
    temp_id.append(c_ids[2])
    ctrl_pts.append(deepcopy(temp))
    ctrl_pts_ids.append(deepcopy(temp_id))

    return ctrl_pts, ctrl_pts_ids


def set_inner_ctrl_pts(facets, layer, edge_pts, num_cpts, change, static):
    """
    Sets the inner control points of the given layer using the edge control points. Structure is similar to
    how the edge pairs were set with different math calculations for determining if the calculate point lies
    within the plane of the facet. Point is determined to be within the triangle boundary using Barycentric coordinates.

    :param facets: List of class type Facet
    :param layer: List of index values associated with facets
    :param edge_pts: List of control point edge pairs
    :param num_cpts: The number of control points to be set between edge pairs
    :param change: Coordinate to be incremented when setting control points (0 - x, 1 - y)
    :param static: Coordinate to left static when setting control points (0 - x, 1 - y)
    :return: A 3-dimensional list with the first dimension being the same length as edge_pts,
            the second dimension being the control points set between the edge pairs, and the
            third dimension simply being the x, y, z values.
            List is constructed this way, in order to more easily construct a 2D grid of all control points later.
    """
    inner_ctrl_pts = []  # Holds the inner control pts
    inner_ctrl_pts_id = []

    for pts in edge_pts:

        # Sets needed variables based on the inputs given
        pos = pts[0][change]
        static_coord = pts[0][static]
        end = pts[1][change]
        spacing = (end - pos) / num_cpts

        inner_pts = []  # Holds the points to be set between the edge pts
        inner_pts_id = []

        # while loops runs the length of the layer in the x-direction
        while pos < end:

            # pos is already incremented as edge pts are already stored
            pos += spacing

            # Used to prevent the distance between final point and the edge from being unnecessarily small
            if (pos + (spacing * 0.2)) > end:
                break

            for f in layer:
                # Calculate the plane made by the points; maintain boundary conditions
                # Determine if coordinates x and y lie within this plane
                # If they are, set z-coordinate to relevant position in relation to x and y
                # Append all 3 coordinates as a list into control_pts and break loop

                # Needed Facet data for calculations, stored into variables for readability and coding
                p0 = facets[f].facet[0]  # Point 0
                p1 = facets[f].facet[1]  # Point 1
                p2 = facets[f].facet[2]  # Point 2
                n = np.cross((p1 - p0), (p2 - p0))  # The normal vector

                # Eq. of a plane is ax + by + cz + d = 0, where (a,b,c) is the normal vector
                d = (n[0] * -p0[0]) + (n[1] * -p0[1]) + (n[2] * -p0[2])

                # Solve for z of the control point
                z = (-(n[change] * pos) - (n[static] * static_coord) - d) / n[2]

                # The control point is constructed first as a list in order to properly build it in [x, y, z] format.
                # Then is converted to an array for mathematical purposes.
                if change == 0:
                    control_pt = np.array([pos, static_coord, z])
                else:
                    control_pt = np.array([static_coord, pos, z])

                # Calculates the UVW coordinates in order to calculate the Barycentric coordinates
                u = p0 - control_pt
                v = p1 - control_pt
                w = p2 - control_pt

                # Calculates the Barycentric coordinates and puts them into an array for mathematical purposes
                alpha = round(np.dot(np.cross(v, w), n) / np.dot(n, n), 3)
                beta = round(np.dot(np.cross(w, u), n) / np.dot(n, n), 3)
                gamma = 1 - alpha - beta
                barycentric_coord = np.array([alpha, beta, gamma])

                # If all the Bary. coords. are between or equal to 0 and 1,
                # then the point is within the facet and is then set
                if (0 <= barycentric_coord).all() and (barycentric_coord <= 1).all():
                    inner_pts.append(control_pt.tolist())
                    inner_pts_id.append(f)
                    break

        inner_ctrl_pts.append(inner_pts)
        inner_ctrl_pts_id.append(inner_pts_id)

    return inner_ctrl_pts, inner_ctrl_pts_id


def inner_pts_exception_case(facets, layer, spacing_x, spacing_y, num_pts_x, num_pts_y, corners):
    """
    Creates all the control points except for the corners. This method is only works for models
    that are a perfect rectangle from the top view.

    :param facets: List of class type Facet
    :param layer: List of index values associated with facets
    :param spacing_x: Spacing of control points in the x-direction
    :param spacing_y: Spacing of control points in the y-direction
    :param num_pts_x: Number of control points in the x-direction
    :param num_pts_y: Number of control points in the y-direction
    :param corners: List of corners in the layer
    :return: Returns a list of the control points
    """

    inner_ctrl_pts = []

    for i in range(num_pts_x + 1):
        x = corners[0][0] + (spacing_x * i)
        y = 0
        for j in range(num_pts_y + 1):
            y = corners[0][1] + (spacing_y * j)
            # a += 1
            if (x + spacing_x) > corners[3][0]:
                x = corners[3][0]
            if (y + spacing_y) > corners[1][1]:
                y = corners[1][1]
            # if x and y in corners:
            # continue
            for f in layer:
                # Calculate the plane made by the points; maintain boundary conditions
                # Determine if coordinates x and y lie within this plane
                # If they are, set z-coordinate to relevant position in relation to x and y
                # Append all 3 coordinates as a list into control_pts and break loop

                # Needed Facet data for calculations, stored into variables for readability and coding
                p0 = facets[f].facet[0]  # Point 0
                p1 = facets[f].facet[1]  # Point 1
                p2 = facets[f].facet[2]  # Point 2
                n = np.cross((p1 - p0), (p2 - p0))  # The normal vector

                # Eq. of a plane is ax + by + cz + d = 0, where (a,b,c) is the normal vector
                d = (n[0] * -p0[0]) + (n[1] * -p0[1]) + (n[2] * -p0[2])

                # Solve for z of the control point
                z = (-(n[0] * x) - (n[1] * y) - d) / n[2]

                control_pt = np.array([x, y, z])

                u = p0 - control_pt
                v = p1 - control_pt
                w = p2 - control_pt

                alpha = round(np.dot(np.cross(v, w), n) / np.dot(n, n), 3)
                beta = round(np.dot(np.cross(w, u), n) / np.dot(n, n), 3)
                gamma = 1 - alpha - beta
                barycentric_coord = np.array([alpha, beta, gamma])

                if (0 <= barycentric_coord).all() and (barycentric_coord <= 1).all():
                    inner_ctrl_pts.append(control_pt)

    return inner_ctrl_pts


def calc_num_pts(corners):
    """
    Calculates the spacing used between control points as if it were a flat plane. Max spacing used is 5 mm.

    :param corners: A sorted list of the corners.
    :return: The spacing of the control points in both the x- and y-direction.
    """

    # Calculating the side lengths using the corner data. All sides are calculated in case of non-equal parallel sides.
    # To prevent a bug that results in a positive result when it should be negative,
    # the absolute value of the difference is taken. This does not change the final result
    diff_x1 = abs(corners[3][0] - corners[0][0])
    diff_x2 = abs(corners[2][0] - corners[1][0])
    diff_y1 = abs(corners[1][1] - corners[0][1])
    diff_y2 = abs(corners[2][1] - corners[3][1])

    # The average of the 2 sides is calculated to keep variables limited, even if there are non-equal parallel sides.
    diff_x_avg = (diff_x1 + diff_x2) / 2
    diff_y_avg = (diff_y1 + diff_y2) / 2

    base_x = round(math.log10(diff_x_avg))
    base_y = round(math.log10(diff_y_avg))

    multiplicity_x = math.pow(10, (base_x - 2))
    multiplicity_y = math.pow(10, (base_y - 2))

    # The total number of points is first divided by 5, as this is the maximum allowed spacing distance.
    # This result is then rounded up to provide a whole number, since a partial point cannot be made.
    num_pts_x = math.ceil(diff_x_avg / (MAX_CTRL_PTS * multiplicity_x))
    num_pts_y = math.ceil(diff_y_avg / (MAX_CTRL_PTS * multiplicity_y))

    # With the number of points now determined, a new spacing distance is calculate using this number.
    spacing_x = diff_x_avg / num_pts_x
    spacing_y = diff_y_avg / num_pts_y

    return num_pts_x, num_pts_y


def evaluate_surface(control_pts):
    """
    Evaluates the NURBS surface using the set control points and global NURBS_DEGREES variable.

    :param control_pts: List of control points to form the surface.
    :return: The NURBS surface after being evaluated and tessellated.
    """

    # Create a BSpline surface instance
    surf = BSpline.Surface()

    # Set degrees
    surf.degree_u = NURBS_DEGREES
    surf.degree_v = NURBS_DEGREES

    # Set control points
    surf.ctrlpts2d = control_pts

    # Set knot vectors
    surf.knotvector_u = utilities.generate_knot_vector(surf.degree_u, surf.ctrlpts_size_u)
    surf.knotvector_v = utilities.generate_knot_vector(surf.degree_v, surf.ctrlpts_size_v)

    # Set evaluation delta
    surf.delta = 0.025

    # Evaluate surface
    surf.evaluate()

    # Sets tessellation algorithm to triangular and then computes the tessellation
    surf.tessellator = tessellate.TriangularTessellate()
    surf.tessellate()

    return surf


def generate_tool_paths(facets, layer, side_surfaces):
    edges = find_edges_id(facets, layer, side_surfaces)
    perimeter = make_perimeter(facets, layer, edges)
    infill = make_infill(perimeter, layer, facets)
    return perimeter + infill


def find_edges_math(facets, layer, side_surfaces):
    edges = []  # The list to hold the line segments that make up the edge
    edge_order = [(0, 1), (0, 2), (1, 2)]

    for f1 in layer:
        found = False
        for e in edge_order:
            edge = [f1.facet[e[0]], f1.facet[e[1]]]
            for index in side_surfaces:
                f2 = facets[index]
                if is_coplanar(edge, f2):
                    edges.append(deepcopy(edge))
                    found = True
                    break
            if found:
                break
    """
    loop: Take facet from layer
        loop through the edges of the facet [(0,1), (0,2), (1,2)]
            loop: check if edge is coplanar to any facet in side_surfaces
                if true:
                    add to line segment to edges
                    break to outer most loop
            
    """

    return edges


def is_coplanar(edge, plane):
    coplanar = True

    # Needed Facet data for calculations, stored into variables for readability and coding
    p0 = plane.facet[0]  # Point 0
    p1 = plane.facet[1]  # Point 1
    p2 = plane.facet[2]  # Point 2
    n = np.cross((p1 - p0), (p2 - p0))  # The normal vector

    # Eq. of a plane is ax + by + cz + d = 0, where (a,b,c) is the normal vector
    d = (n[0] * -p0[0]) + (n[1] * -p0[1]) + (n[2] * -p0[2])

    for pt in edge:
        store = (n[0] * pt[0]) + (n[0] * pt[1]) + (n[0] * pt[2]) + d
        if (n[0] * pt[0] + n[0] * pt[1] + n[0] * pt[2] + d) != 0:
            coplanar = False
            # break

    return coplanar


def find_edges_id(facets, layer, side_surfaces):
    """
    Finds the edges of the layer.

    :param facets: A list of all the facets and relevant data.
    :param layer: The list of indexes that make up the layer.
    :param side_surfaces: The list of indexes that make up the side surfaces.
    :return: A list of the line segments that make up the edge.
    """

    edges = []  # The list to hold the line segments that make up the edge

    # Gets the index in layer
    for f1 in layer:

        # Using the index, get the key to the adjacent dictionary associated with that facet
        for key in facets[f1].adjacent:

            # Gets the index in side_surfaces
            for f2 in side_surfaces:

                # Check if the adjacent dictionary value is the same as the index value from side_surfaces
                if facets[f1].adjacent[key] == f2:
                    # Adds the side data to edges list via the use of the key
                    edge = Edge(facets[f1].sides[key], f1)
                    edges.append(edge)
                    break
    return edges


def make_perimeter(facets, layer, edges):
    perimeter = []
    for edge in edges:
        pass
    return perimeter


def make_infill(perimeter, layer, facets):
    infill = []

    return infill


# This provided line is required at the end of a Python file
# to call the main() function.
if __name__ == '__main__':
    main()
