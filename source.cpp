#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "library.h"

using namespace std;

// Stores information about a named place in the United States.
struct city {
	string state, name;
	double latitude, longitude;
	int intersect;
	city * next;

	city(string s, string n, double la, double lo, int i) {
		state = s; name = n; latitude = la; longitude = lo;
		intersect = i; next = NULL;
	}
};

// Turns a line from named-places.txt into a city struct with the relevant data.
city * process_line(string line) {
	string state = line.substr(8, 2);
	string name = line.substr(10, line.find("  ") - 10);
    double latitude = atof(line.substr(80, 10).c_str());
	double longitude = atof(line.substr(90, 11).c_str());
	int intersect = atoi(line.substr(101, 5).c_str());
	return new city(state, name, latitude, longitude, intersect);
}

// Hashes a string.
int hash_string(string s, int hsize) {
    int h = 7761;
    for (int i = 0; i < s.length(); i++)
        h *= (s[i] + 7);
    h = h % hsize;
    if (h < 0)
        h = -h;
    return h;
}

// Reads in cities from named-places.txt and populates a hashtable with them.
void populate_hashtable(city ** hashtable) {
    for (int i = 0; i < 5000; i++)
        hashtable[i] = NULL;

    ifstream file("named-places.txt");
	if (!file.is_open()) {
		cout << "File not found: named-places.txt.\n";
		return;
	}

    string line;
    while (getline(file, line)) {
		// Most of the heavy lifting is done by process_line().
        city * newplace = process_line(line);
        int newhash = hash_string(newplace->name, 5000);
		// I used closed addressing, so each city struct will either be added directly
		// to the hashtable, or it will be appended to the end of the last city struct
		// that is already present.
        if (hashtable[newhash] == NULL)
            hashtable[newhash] = newplace;
        else {
            city * current = hashtable[newhash];
            while (current->next != NULL)
                current = current->next;
            current->next = newplace;
        }
    }
}

struct intersection;
// A road connecting two intersections.
struct connection {
	string name;
	intersection * A, * B;
	double length;

	connection(string n, intersection * first, intersection * last, double l) {
		name = n;
		A = first;
		B = last;
		length = l;
	}

	// Passed a pointer to an intersection, will return the other intersection
	// it is connected to.
	intersection * other(intersection * cur) {
		if (A == cur)
			return B;
		return A;
	}
};

// A place where two roads meet.
struct intersection {
	string state, city;
	int num;
	double dist, shortest, latitude, longitude;
	vector<connection*> connections;

	intersection(string s, string c, int n, double di, double la, double lo) {
		state = s;
		city = c;
		num = n;
		dist = di;
		shortest = -1.0;
		latitude = la;
		longitude = lo;
	}

	// Passed an intersection, will return the road connecting the two.
	connection * getconnection(intersection * other) {
		for (int i = 0; i < connections.size(); i++)
			if (connections[i]->A == other || connections[i]->B == other)
				return connections[i];
		return NULL;
	}
};

// Used for running the priority queue.
struct node {
    intersection * s;
    double p;
    node(intersection * st, double pr) {
        s = st;
        p = pr;
    }
};

// A simple priority queue. Has 3 functions: add(), removel() and empty().
class priorityqueue {
private:
    vector<node*> arr;
public:
	bool empty() {
		return arr.size() == 0;
	}

    void add(intersection * s, double prio) {
        arr.push_back(new node(s, prio));
		// Iterate up the tree, swapping the new item w/ parent if it has a higher priority.
        if (arr.size() > 1) {
            int i = arr.size() - 1;
            for (int p = i / 2; p >= 0; p /= 2) {
                if (arr[p]->p > arr[i]->p)
                    swap(arr[p], arr[i]);
                else
                    break;
                i /= 2;
                if (p == 0)
                    break;
            }
        }
    }

	// Remove the item with the highest priority.
    intersection * removel() {
        int i = arr.size() - 1;
		if (i == -1) // Make sure it's not empty.
			return NULL;
		else if (i == 0) { // If there's one item, return it and we're done.
			intersection * s = arr[i]->s;
			arr.pop_back();
			return s;
		}
		// If there's two items or more, swap the highest priority item with the
		// last item in the heap. Pop back the highest priority item to be returned
		// later.
        swap(arr[0], arr[i]);
        intersection * s = arr[i]->s;
        arr.pop_back();

		i--;

        int curr = 0, p1, p2;
        while (true) {
			// Since the lowest item has been moved to the front, it's probably out of
			// place. It is checked against the children if it needs to be swapped.
			p1 = curr * 2 + 1;
			p2 = curr * 2 + 2;

			// If there aren't any children, then the item is fine where it is.
			if (curr >= i || p1 >= i)
				break;
			// If there's one child and the child has a higher priority, swap it.
			else if (p2 >= i && arr[curr]->p > arr[p1]->p) {
				swap(arr[curr], arr[p1]);
				break;
			} else {
				// If there's two children, find the lowest. If it's lower than the current,
				// swap them. If not, break.
				if (arr[p1]->p >= arr[p2]->p && arr[curr]->p > arr[p2]->p) {
					swap(arr[curr], arr[p2]);
					curr = p2;
				} else if (arr[p2]->p > arr[p1]->p && arr[curr]->p > arr[p1]->p) {
					swap(arr[curr], arr[p1]);
					curr = p1;
				} else {
					break;
				}
			}
        }
        return s;
    }
};

// Reads intersections from intersections.txt and adds them to the array passed.
void read_intersections(intersection ** intersections) {
	ifstream file("intersections.txt", ios::in);
	if (!file.is_open()) {
		cout << "Fatal error: Couldn't open intersections.txt.\n";
		return;
	}
	string line;
	int index = 0;
	// Read in the intersection data and add it to intersections.
    while (getline(file, line)) {
        double longitude = atof(line.substr(0, 9).c_str());
        double latitude = atof(line.substr(12, 7).c_str());
        double dist = atof(line.substr(21, 6).c_str());
        string state = line.substr(28, 2);
        string city = line.substr(31, line.size()-31);
        intersections[index] = new intersection(state, city, index, dist, latitude, longitude);
        index++;
    }
}

// Reads in connections from connections.txt.
void read_connections(intersection ** intersections) {
	string line;
    ifstream file("connections.txt", ios::in);
	if (!file.is_open()) {
		cout << "Fatal error: Couldn't open connections.txt.\n";
		return;
	}
    while (getline(file, line)) {
		// I used a stringstream because I felt it was less ugly than calling string.find(" ") a
		// bajillion times.
        stringstream iss(line);
        string part, name;
        int A, B, counter = 0;
        double dist;
        while (getline(iss, part, ' ')) {
            if (counter == 0)
                name = part;
            else if (counter == 2)
                A = atoi(part.c_str());
            else if (counter == 3)
                B = atoi(part.c_str());
            else if (counter == 4)
                dist = atof(part.c_str());
            counter++;
        }

		// Create a new connection and add it to the intersection structs it connects.
        connection * c = new connection(name, intersections[A], intersections[B], dist);
        intersections[A]->connections.push_back(c);
        intersections[B]->connections.push_back(c);
    }
}

// Gets the user input for two cities. These cities will be navigated between.
void read_in(string * c1, string * s1, string * c2, string * s2) {
	cout << "Enter start city: ";
	getline(cin, *c1);
	cout << "Enter start state (abbr): ";
	getline(cin, *s1);
	cout << "Enter destination city: ";
	getline(cin, *c2);
	cout << "Enter destination state (abbr): ";
	getline(cin, *s2);
}

// Is passed the cities the user input in read_in(). Looks them up in
// the hashtable to find the nearest intersection. Returns the number of
// the intersection so the shortest path function can do its job.
int lookup_intersection(city ** hashtable, string c, string s) {
	int h = hash_string(c, 5000);
	city * cur = hashtable[h];
	while (cur != NULL) {
		if (cur->name == c && cur->state == s)
			return cur->intersect;
		cur = cur->next;
	}
	return -1;
}

// Finds the shortest path between intersection #'s start and finish. Intersections is the array containing
// all of the intersections. Route is an empty vector which will be filled and used later to print directions and
// draw them on the map.
void find_shortest(intersection ** intersections, int start, int finish, vector<intersection*> &route) {
	priorityqueue pq;
	intersections[start]->shortest = 0.0;
    pq.add(intersections[start], 0.0); // Add the start and get the priority queue going.

    while (!pq.empty()) {
        intersection * curr = pq.removel(); // Remove the lowest from the priority queue.
		if (curr == NULL)
			break;
		// Go through each of the connections/roads going out of it.
        for (int i = 0; i < curr->connections.size(); i++) {
            connection * r = curr->connections[i];
            double newlength = curr->shortest + r->length;
			// Find the distance to the intersection at the other end. If the intersection at the other
			// end hasn't been checked yet (shortest = -1) or has a shortest distance higher than the one
			// we've just found, add it to the priority queue and adjust its shortest.

			// Also, quit searching on this path if we've already found a shorter path to the finish.
			// If there's already a 25 mile path to the end, there's no need to continue if we're starting
			// at 30 miles from another place.
            if (intersections[finish]->shortest != -1 && newlength > intersections[finish]->shortest)
                continue;
			intersection * endpoint = r->other(curr);
            if (endpoint->shortest == -1 || (endpoint->shortest != -1 && endpoint->shortest > newlength)) {
                endpoint->shortest = newlength;
				// Only add it to the priority queue if it's not the final point.
                if (endpoint != intersections[finish])
                    pq.add(endpoint, newlength);
			}
        }
    }

	// If the finish has a shortest value of -1, that means the intersections aren't connected.
    if (intersections[finish]->shortest == -1)
		return;

	// Now work through the intersections backwards to find the route we took to get the shortest path.
	intersection * curr = intersections[finish];
    while (curr != intersections[start]) {
        for (int i = 0; i < curr->connections.size(); i++) {
            connection * r = curr->connections[i];
            intersection * endpoint = NULL;
			endpoint = r->other(curr);
			// Since we're working with doubles they probably won't match up exactly. As long
			// as they're within 0.1 of each other it should be the right node.
            double x = endpoint->shortest - (curr->shortest - r->length);
            if (x < 0.1 && x > -0.1) {
                route.push_back(curr);
                curr = endpoint;
                break;
            }
        }
    }
	// The route vector is finish -> start so it needs to be reversed.
    route.push_back(intersections[start]);
	reverse(route.begin(), route.end());
	// No need to return anything since we're altering the original route vector.
}

// Returns the direction traveling from A to B (e.g. North, Southeast).
string direction(intersection * A, intersection * B) {
	double lat = B->latitude - A->latitude;
	double lon = B->longitude - A->longitude;
	double theta = 180.0 * atan(lat / lon) / 3.14159265;
	
    if (lat >= 0 && lon >= 0) // Quadrant 1
        if (theta < 22.5)
            return "East";
        else if (theta < 67.5)
            return "Northeast";
        else
            return "North";
    else if (lat >= 0 && lon <= 0) // Quadrant 2
        if (-theta < 22.5)
            return "North";
        else if (-theta < 67.5)
			return "Northwest";
        else
            return "West";
    else if (lat <= 0 && lon <= 0) // Quadrant 3
        if (theta < 22.5)
            return "West";
        else if (theta < 67.5)
            return "Southwest";
        else
            return "South";
    else if (lat <= 0 && lon >= 0) // Quadrant 4
        if (-theta < 22.5)
            return "South";
        else if (-theta < 67.5)
            return "Southeast";
        else
            return "East";
	return "";
}

// Passed the route obtained in shortest_path(), prints out directions. Will only
// print when the road name changes.
void print_directions(vector<intersection*> &route) {
	double sum = 0.0; // Keep a tally of how far we've gone on the current road.
	for (int i = 0; i < route.size() - 1; i++) {
		// Only need to iterate through the end minus one, as there's nowhere to go
		// once the end is reached.
		string dir = direction(route[i], route[i+1]);
		connection * c = route[i]->getconnection(route[i+1]);
		connection * c2 = NULL;
		if (i+2 < route.size())
			c2 = route[i+1]->getconnection(route[i+2]);
		sum += c->length;
		bool is_same = false;
		if (c2 != NULL)
			if (c->name.find(c2->name) != string::npos || c2->name.find(c->name) != string::npos)
				is_same = true;
		if (is_same) {
			continue;
		} else {
			cout << dir << " on " << c->name << " for " << sum << " miles to near ";
			cout << route[i+1]->city << ", " << route[i+1]->state << "\n";
			sum = 0.0;
		}
	}
}

// Stores information about a map.dat from coverage.txt.
struct mapfile {
	int lat_north, lat_south, lon_west, lon_east, d;
	string filename;
	mapfile(int ln, int ls, int lw, int le, string f, int d2) {
		lat_north = ln;
		lat_south = ls;
		lon_west = lw;
		lon_east = le;
		filename = f;
		d = d2;
	}
};

// Passed the route we obtained in shortest_path(), gets the minimal map size needed
// to display the route.
mapfile * get_map(vector<intersection*> &route) {
	// Find the largest and smallest latitudes and longitudes.
	double lat_smallest = route[0]->latitude,
		   lat_largest = lat_smallest,
		   lon_smallest = route[0]->longitude,
		   lon_largest = lon_smallest;
	for (int i = 1; i < route.size(); i++) {
		double lat = route[i]->latitude;
		double lon = route[i]->longitude;
		if (lat < lat_smallest)
			lat_smallest = lat;
		if (lat > lat_largest)
			lat_largest = lat;
		if (lon < lon_smallest)
			lon_smallest = lon;
		if (lon > lon_largest)
			lon_largest = lon;
	}

	// Read coverage.txt and create an array storing all of the info about the map files.
	mapfile * mapfiles[149];
	ifstream file("coverage.txt", ios::in);
	if (!file.is_open()) {
		cout << "Fatal error: Couldn't find coverage.txt.\n";
		return NULL;
	}

	string line;
	int index = 0;
	while (getline(file, line)) {
		int lat_north, lat_south, lon_west, lon_east, c = 0;
		// I use stringstream here also because I don't feel like calling
		// string.find() a bunch of times.
		string filename, s;
		stringstream iss(line);
		while (getline(iss, s, ' ')) {
			if (c == 0)
				lat_north = atof(s.c_str());
			else if (c == 1)
				lat_south = atof(s.c_str());
			else if (c == 2)
				lon_west = atof(s.c_str());
			else if (c == 3)
				lon_east = atof(s.c_str());
			else
				filename = s;
			c++;
		}
		// Get the size of the map, D5 -> 5, D30 -> 30, etc.
		string d = filename.substr(filename.length()-6, 2);
		if (d.substr(0, 1) == "D")
			d = d.substr(1);
		int di = atof(d.c_str());

		mapfiles[index] = new mapfile(lat_north, lat_south, lon_west, lon_east, filename, di);
		index++;
	}

	// Locate smallest tile containing the necessary latitudes and longitudes. Does a pass
	// for D = 5, and then every multiple of 10 up to 80. Breaks as soon as it finds a map
	// that works, so it will return the smallest necessary map.
	int d = 5;
	mapfile * themap = NULL;
	// If we hit 80 and we haven't found anything (which shouldn't happen), break so we're
	// not stuck in an infinite loop.
	while (d < 80) {
		// There are 148 map files so check each one.
		for (int i = 0; i < 149; i++) {
			mapfile * m = mapfiles[i];
			// Only check them if they're the correct size for what pass we're on.
			if (m->d != d)
				continue;
			// If all of the latitudes and longitudes will fit, we're done. Return the map.
			if (lat_smallest >= m->lat_south && lat_largest <= m->lat_north &&
				lon_smallest >= m->lon_west && lon_largest <= m->lon_east) {
				return m;
			}
		}
		// As above, if d is 5, it becomes 10, otherwise just add 10. I don't think there are any
		// maps of size 40, 50, or 70, but it's simpler to do it this way.
		if (d == 5)
			d = 10;
		else
			d += 10;
	}

	if (themap == NULL) {
		cout << "Couldn't find a map that fit the points.\n";
	}
	return themap;
}


// Draws the map and the route.
void draw_map(mapfile * map, vector<intersection*> & route) {
	ifstream file("map/" + map->filename, ios::in|ios::binary);
	if (!file.is_open()) {
		cout << "Fatal error: Couldn't find " << map << ".\n";
		return;
	}

	// It took me a bit to realize the maps weren't all square. Adjust the size of the window
	// and how we draw based on the map size.
	int rmax, cmax;
	if (map->d <= 30) {
		rmax = 600;
		cmax = 600;
	} else if (map->d == 60) {
		rmax = 400;
		cmax = 800;
	} else {
		rmax = 500;
		cmax = 800;
	}

	short int * s = new short int;
	int row = 0, column = 0;
	// Move it to (50, 50) for convenience.
	make_window(cmax, rmax, 50, 50);

	file.read((char*)s, 2);
	while (!file.eof()) {
		file.read((char*)s, 2);
		int c;

		if (*s == -500)
			c = color::cyan;
		else {
			// I found that using 600.0 as a scale made a map that looked nice.
			float b = * s / 600.0;
			c = color::rgb(0, b, 0);
		}
		set_pixel_color(column-1, row, c);

		// Move down the column and then down the row once you reach the end.
		column += 1;
		if (column > cmax-1) {
			column = 0;
			row += 1;
		}
	}

	// Draw the route.
	set_pen_width(5);
	set_pen_color(color::red);

	for (int i = 0; i < route.size() - 1; i++) {
		/*
		Points are translated from (longitude, latitude) to (x, y) that fits on the
		drawing window. To do this: subtract the west longitude from x and the
		north latitude from y so that the upper corners line up. Then, the point has
		to be scaled so the rest of the points will fit on the map. For x the scale
		is:

		- (height / (longitude_west - longitude_east))
		*/
		double x1 = route[i]->longitude, x2 = route[i+1]->longitude,
			y1 = route[i]->latitude, y2 = route[i+1]->latitude;
		x1 = (x1 - map->lon_west) * -(cmax / (map->lon_west - map->lon_east));
		x2 = (x2 - map->lon_west) * -(cmax / (map->lon_west - map->lon_east));
		y1 = (y1 - map->lat_north) * -(rmax / (map->lat_north - map->lat_south));
		y2 = (y2 - map->lat_north) * -(rmax / (map->lat_north - map->lat_south));

		// Draw from point A to point B to illustrate the route.
		move_to(x1, y1);
		draw_to(x2, y2);
	}
}

void main() {
	// Create the hashtable of cities.
    city * hashtable[5000];
	for (int i = 0; i < 5000; i++)
		hashtable[i] = NULL;
	populate_hashtable(hashtable);

	// Get cities to navigate from user and locate nearest intersections.
	string c1, s1, c2, s2;
	read_in(&c1, &s1, &c2, &s2);
	int start = lookup_intersection(hashtable, c1, s1), finish = lookup_intersection(hashtable, c2, s2);
	cout << start << " " << finish << "\n";
	if (start == -1 && finish == -1) {
		cout << "Neither of the places you've entered seems to exist. Exiting.\n";
		return;
	} else if (start == -1) {
		cout << "Your starting place doesn't seem to exist. Exiting.\n";
		return;
	} else if (finish == -1) {
		cout << "Your ending place doesn't seem to exist. Exiting.\n";
		return;
	}

	// Read in intersections and roads.
	intersection * intersections[29146];
	read_intersections(intersections);
	read_connections(intersections);

	// Find shortest path.
	vector<intersection*> route;
	find_shortest(intersections, start, finish, route);

	// Print directions.
	print_directions(route);

	// Draw the map.
	mapfile * map = get_map(route);
	draw_map(map, route);
}