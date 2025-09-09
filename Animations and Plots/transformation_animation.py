from manim import *

from TerraFrame import CelestialTerrestrial
from TerraFrame.Utilities import Conversions
from TerraFrame.Utilities import TransformationMatrices
from TerraFrame.Utilities.Time import JulianDate


def get_matrices():
    jd_utc = JulianDate.JulianDate.j2000() + 500.7
    jd_tt = Conversions.utc_to_tt(jd_utc)

    ct = CelestialTerrestrial.CelestialTerrestrialTransformation()

    ct.gcrs_to_itrs(jd_tt)

    return ct.t_gi, ct.t_gc, ct.t_ct, ct.t_ti


class EarthFrameAnimation(ThreeDScene):
    def construct(self):
        # Set up 3D camera
        self.set_camera_orientation(phi=70 * DEGREES, theta=30 * DEGREES)

        # Baseline reference axes
        base_axes = ThreeDAxes(x_range=[-1, 1, 1], y_range=[-1, 1, 1],
                               z_range=[-1, 1, 1], axis_config={"color": WHITE})

        # Second (rotating) axes
        rot_axes = base_axes.copy()

        x_axis_label = Text('X')
        # noinspection PyUnresolvedReferences
        x_axis_label.next_to(base_axes.x_axis.get_end(), RIGHT)
        x_axis_label.rotate(90 * DEGREES, axis=X_AXIS)
        x_axis_label.rotate(140 * DEGREES, axis=Z_AXIS)
        x_axis_label.rotate(-15 * DEGREES, axis=OUT)

        y_axis_label = Text('Y')
        # noinspection PyUnresolvedReferences
        y_axis_label.next_to(base_axes.y_axis.get_end(), UP)
        y_axis_label.rotate(90 * DEGREES, axis=X_AXIS)
        y_axis_label.rotate(100 * DEGREES, axis=Z_AXIS)
        y_axis_label.rotate(20 * DEGREES, axis=OUT)

        z_axis_label = Text('Z')
        # noinspection PyUnresolvedReferences
        z_axis_label.next_to(base_axes.z_axis.get_end(), UP)
        z_axis_label.rotate(90 * DEGREES, axis=X_AXIS)
        z_axis_label.rotate(135 * DEGREES, axis=Z_AXIS)
        z_axis_label.rotate(-10 * DEGREES, axis=OUT)

        self.add(base_axes, x_axis_label, y_axis_label, z_axis_label)

        # Sphere fixed to second axes
        sphere = Sphere(radius=1.5, resolution=(24, 48), fill_opacity=1.0,
                        checkerboard_colors=[GREEN_D, GREEN_E])

        equator = Circle(radius=1.5, color=BLACK)

        group = VGroup(rot_axes, sphere, equator)

        # Title setup
        title = VGroup()
        t0 = Text("CGRS", color=WHITE).scale(0.5)
        title.add(t0)
        self.add_fixed_in_frame_mobjects(title)
        title.to_corner(UL)

        date_note = Text("2001-05-16 04:48:00 UTC", color=WHITE)
        date_note.scale(0.5)
        self.add_fixed_in_frame_mobjects(date_note)
        date_note.to_corner(DL)

        def extend_title(label, text_color):
            arrow = Text(" ⟶ ").scale(0.5).next_to(title[-1], RIGHT, buff=0.1)
            word = Text(label, color=text_color).scale(0.5).next_to(arrow,
                                                                    RIGHT,
                                                                    buff=0.1)
            title.add(arrow, word)
            self.add_fixed_in_frame_mobjects(title)
            title.to_corner(UL)
            return [FadeIn(arrow), FadeIn(word)]

        self.play(Create(base_axes), Create(rot_axes), Create(x_axis_label),
                  Create(y_axis_label), Create(z_axis_label), Create(title),
                  Create(date_note))
        self.play(Create(sphere), Create(equator))
        self.wait(2)

        _, t_gc, t_ct, t_ti = get_matrices()

        angle, axis = (
            TransformationMatrices.angle_and_axis_from_transformation(t_gc.T))

        note = Text("Angle scaled 1,000x", color=RED)
        note.scale(0.5)
        self.add_fixed_in_frame_mobjects(note)
        note.to_corner(DR)

        # Perform rotation and color fade simultaneously
        new_text_animated = extend_title("CIRS", RED)
        self.play(rot_axes.animate.set_color(RED).set_run_time(1),
                  *new_text_animated, Create(note))
        self.play(Rotate(group, angle=angle * 1e3, axis=axis, run_time=2))
        self.wait(1)

        rot_old = rot_axes.copy()
        self.add(rot_old)
        angle, axis = (
            TransformationMatrices.angle_and_axis_from_transformation(t_ct.T))
        new_text_animated = extend_title("TIRS", BLUE)
        self.play(rot_axes.animate.set_color(BLUE).set_run_time(1),
                  *new_text_animated, Uncreate(note))
        self.play(Rotate(group, angle=angle, axis=axis), run_time=2)
        self.wait(1)

        note = Text("Angle scaled 100,000x", color=ORANGE)
        note.scale(0.5)
        self.add_fixed_in_frame_mobjects(note)
        note.to_corner(DR)

        rot_old = rot_axes.copy()
        self.add(rot_old)
        angle, axis = (
            TransformationMatrices.angle_and_axis_from_transformation(t_ti.T))
        new_text_animated = extend_title("ITRS", ORANGE)
        self.play(rot_axes.animate.set_color(ORANGE).set_run_time(1),
                  *new_text_animated, Create(note))
        self.play(Rotate(group, angle=angle * 1e5, axis=axis), run_time=2)
        self.wait(3)
