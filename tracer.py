import math
import colorsys
from PIL import Image

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 1000

pixel_values = []

screen = Image.new("RGB", (SCREEN_WIDTH, SCREEN_HEIGHT))


class Vector:
    def __init__(self, x, y, z):
        self.vector = (x, y, z)

    def __str__(self):
        return str(self.vector)

    @property
    def getVector(self):
        return Vector(x, y, z)

    def setVector(self, x, y, z):
        self.vector = (x, y, z)

    def get_x(self):
        return self.vector[0]

    def get_y(self):
        return self.vector[1]

    def get_z(self):
        return self.vector[2]

    def minusVector(self, vector):
        return Vector(self.vector[0] - vector.vector[0], self.vector[1] - vector.vector[1],
                      self.vector[2] - vector.vector[2])

    def addVector(self, vector):
        return Vector(self.vector[0] + vector.vector[0], self.vector[1] + vector.vector[1],
                      self.vector[2] + vector.vector[2])

    def multiplyByScalar(self, scalar):
        return Vector(self.vector[0] * scalar, self.vector[1] * scalar, self.vector[2] - scalar)

    def magnitude(self):
        return math.sqrt(self.vector[0] ** 2 + self.vector[1] ** 2 + self.vector[2] ** 2)

    def normalise(self):
        magnitude = self.magnitude()
        return Vector(self.vector[0] / magnitude, self.vector[1] / magnitude, self.vector[2] / magnitude)

    def dot(self, vector):
        return self.get_x() * vector.get_x() + self.get_y() * vector.get_y() + self.get_z() * vector.get_z()


class Sphere:
    def __init__(self, centre_position, r, colour, phong):
        self.center_position = centre_position
        self.radius = r
        self.colour = colour
        self.phong = phong

    def get_position(self):
        return self.center_position

    def get_radius(self):
        return self.radius

    def get_colour(self):
        return self.colour

    def get_phong_exp(self):
        return self.phong

    def surface_normal_at_point(self, point):
        return point.minusVector(self.center_position).normalise()


class Plane:
    def __init__(self, position, normal, colour, phong):  # r.n = p
        self.position = position
        self.normal = normal
        self.colour = colour
        self.phong = phong

    def get_position(self):
        return self.position

    def get_normal(self):
        return self.normal

    def get_colour(self):
        return self.colour

    def get_phong_exp(self):
        return self.phong

    def get_plane_scalar(self):
        return self.position.get_x() * self.normal.get_x() + \
               self.position.get_y() * self.normal.get_y() + \
               self.position.get_z() * self.normal.get_z()

    def surface_normal_at_point(self, point):
        return self.normal


class Ray:
    def __init__(self, position, direction):
        self.position = position
        self.direction = direction

    def find_sphere_intersection(self, sphere):
        i = self.position.get_x()
        j = self.position.get_y()
        k = self.position.get_z()
        a = self.direction.get_x()
        b = self.direction.get_y()
        c = self.direction.get_z()
        r = sphere.get_radius()
        d = sphere.get_position().get_x()
        e = sphere.get_position().get_y()
        f = sphere.get_position().get_z()
        A = a ** 2 + b ** 2 + c ** 2
        B = 2 * (a * (i - d) + b * (j - e) + c * (k - f))
        C = (i ** 2 + j ** 2 + k ** 2) + (d ** 2 + e ** 2 + f ** 2) - 2 * (i * d + j * e + k * f) - r ** 2
        if num_of_sphere_intersections(A, B, C) == 0:
            return None
        else:
            distance = min(solve_quadratic(A, B, C))
            intersection = Vector(i + distance * a, j + distance * b, k + distance * c)
            return intersection

    def find_plane_intersection(self, plane):
        scalar = plane.get_plane_scalar()
        i = self.position.get_x()
        j = self.position.get_y()
        k = self.position.get_z()
        x = self.direction.get_x()
        y = self.direction.get_y()
        z = self.direction.get_z()
        d = plane.get_normal().get_x()
        e = plane.get_normal().get_y()
        f = plane.get_normal().get_z()
        numerator = scalar - (d * i + e * j + f * k)
        denominator = d * x + e * y + f * z
        try:
            distance = numerator / denominator
        except ZeroDivisionError:
            return None
        if distance > 0:
            intersection = Vector(i + distance * x, j + distance * y, k + distance * z)
            return intersection
        else:
            return None


class Scene:
    def __init__(self):
        self.sceneArr = []
        self.camera_position = Vector(0, 0, 0)
        self.light_sources = []

    def addSphere(self, sphere):
        self.sceneArr.append(sphere)

    def addPlane(self, plane):
        self.sceneArr.append(plane)

    def add_light_source(self, position):
        self.light_sources.append(position)

    def get_light_sources(self):
        return self.light_sources

    def set_camera_position(self, x, y, z):
        self.camera_position = Vector(x, y, z)

    def get_camera_position(self):
        return self.camera_position

    def get_sceneArr(self):
        return self.sceneArr


def solve_quadratic(a, b, c):
    return [((-1 * b) + math.sqrt(b ** 2 - (4 * a * c))) / 2 * a, ((-1 * b) - math.sqrt(b ** 2 - (4 * a * c))) / 2 * a]


def num_of_sphere_intersections(a, b, c):
    discriminant = b ** 2 - (4 * a * c)
    if discriminant == 0:
        return 1
    if discriminant > 0:
        return 2
    if discriminant < 0:
        return 0


def cartesian_to_pixel_Coordinates(x, y):
    screen_x = x + SCREEN_WIDTH / 2
    screen_y = SCREEN_HEIGHT / 2 - y
    return screen_x, screen_y


def pixel_to_cartesian_Coordinates(x, y):
    cart_x = x - SCREEN_WIDTH / 2
    cart_y = SCREEN_HEIGHT / 2 - y
    return cart_x, cart_y


scene = Scene()
FOV = 60
camera_distance = (SCREEN_WIDTH / 2) / math.tan(math.radians(FOV / 2))
scene.set_camera_position(0, 0, -1 * camera_distance)
camera_position = scene.get_camera_position()

adj_width = int(0.05 * SCREEN_WIDTH)
adj_height = int(0.05 * SCREEN_HEIGHT)

scene.add_light_source(Vector(-1 * (SCREEN_WIDTH / 2 - adj_width), SCREEN_HEIGHT / 2 - adj_height,
                              SCREEN_HEIGHT - adj_height))  # back left
# scene.add_light_source(Vector(SCREEN_WIDTH/2 - adj_width, SCREEN_HEIGHT/2 - adj_height, SCREEN_HEIGHT - adj_height))  # back right
# scene.add_light_source(Vector(-1*(SCREEN_WIDTH/2 - adj_width), SCREEN_HEIGHT/2 - adj_height, adj_height))  # front left
scene.add_light_source(Vector(SCREEN_WIDTH/2 - adj_width, SCREEN_HEIGHT/2 - adj_height, adj_height))  # front right
light_sources = scene.get_light_sources()

sphere1 = Sphere(Vector(-300, -300, 200), 200, (200, 40, 155), 100)
scene.addSphere(sphere1)
sphere2 = Sphere(Vector(300, -450, 100), 50, (20, 40, 205), 100)
scene.addSphere(sphere2)

plane1 = Plane(Vector(0, -1 * SCREEN_HEIGHT / 2, 0), Vector(0, 1, 0).normalise(), (255, 0, 0), 100)  # floor
scene.addPlane(plane1)
plane2 = Plane(Vector(0, 0, SCREEN_HEIGHT), Vector(0, 0, -1).normalise(), (0, 255, 0), 100)  # back
scene.addPlane(plane2)
plane3 = Plane(Vector(-1 * SCREEN_WIDTH / 2, 0, 0), Vector(1, 0, 0).normalise(), (0, 0, 255), 100)  # left
scene.addPlane(plane3)
plane4 = Plane(Vector(SCREEN_WIDTH / 2, 0, 0), Vector(-1, 0, 0).normalise(), (0, 0, 255), 100)  # right
scene.addPlane(plane4)
plane5 = Plane(Vector(0, SCREEN_HEIGHT / 2, 0), Vector(0, -1, 0).normalise(), (15, 236, 252), 100)  # ceiling
scene.addPlane(plane5)

intensity = 1


def shadow_ray_intersects(intersection_point, shadow_ray_dir, distance):
    for scene_object in scene.get_sceneArr():
        shadow_ray = Ray(intersection_point, shadow_ray_dir)
        if type(scene_object) is Sphere:
            intersection = shadow_ray.find_sphere_intersection(scene_object)
            if intersection is None:
                continue
            else:
                closeness = intersection.minusVector(intersection_point).magnitude()
                if closeness < distance:
                    return True
        if type(scene_object) is Plane:
            intersection = shadow_ray.find_plane_intersection(scene_object)
            if intersection is None:
                continue
            else:
                closeness = intersection.minusVector(intersection_point).magnitude()
                if closeness < distance:
                    return True


for i in range(SCREEN_HEIGHT):
    for j in range(SCREEN_WIDTH):
        cart_coordinates = pixel_to_cartesian_Coordinates(j, i)
        screen_coordinates = Vector(cart_coordinates[0], cart_coordinates[1], 0)
        ray = Ray(camera_position, screen_coordinates.minusVector(camera_position).normalise())
        closest_object = None
        for item in scene.get_sceneArr():
            if type(item) is Sphere:
                intersection_point = ray.find_sphere_intersection(item)
            elif type(item) is Plane:
                intersection_point = ray.find_plane_intersection(item)
            if intersection_point is not None and closest_object is not None:
                if type(closest_object) is Sphere:
                    closest_intersection = ray.find_sphere_intersection(closest_object)
                elif type(closest_object) is Plane:
                    closest_intersection = ray.find_plane_intersection(closest_object)
                if intersection_point.minusVector(camera_position).magnitude() < closest_intersection.minusVector(
                        camera_position).magnitude():
                    closest_object = item
            elif intersection_point is not None and closest_object is None:
                closest_object = item

        if closest_object is None:
            pixel_values.append((0, 0, 0))
            continue
        item_base_colour = closest_object.get_colour()
        phong_exp = closest_object.get_phong_exp()
        red = 0.2 * item_base_colour[0] / 255
        green = 0.2 * item_base_colour[1] / 255
        blue = 0.2 * item_base_colour[2] / 255
        if type(closest_object) is Sphere:
            intersection_point = ray.find_sphere_intersection(closest_object)
        elif type(closest_object) is Plane:
            intersection_point = ray.find_plane_intersection(closest_object)
        for k in range(len(light_sources)):
            normal = closest_object.surface_normal_at_point(intersection_point).normalise()
            light_dir = light_sources[k].minusVector(intersection_point).normalise()
            if shadow_ray_intersects(intersection_point.addVector(light_dir.multiplyByScalar(0.1)), light_dir,
                                     light_sources[k].minusVector(intersection_point).magnitude()):
                continue
            n_dot_l = normal.dot(light_dir)
            v = camera_position.minusVector(intersection_point).normalise()
            half_vector = v.addVector(light_dir).normalise()
            n_dot_h = normal.dot(half_vector)
            red += intensity * (item_base_colour[0] / 255 * max(0, n_dot_l)) + (
                    intensity * 160 / 255 * (max(0, n_dot_h) ** phong_exp))
            green += intensity * (item_base_colour[1] / 255 * max(0, n_dot_l)) + (
                    intensity * 160 / 255 * (max(0, n_dot_h) ** phong_exp))
            blue += intensity * (item_base_colour[2] / 255 * max(0, n_dot_l)) + (
                    intensity * 160 / 255 * (max(0, n_dot_h) ** phong_exp))
        pixel_values.append((int(red * 255), int(green * 255), int(blue * 255)))

screen.putdata(pixel_values)
screen.save('test.png')
