#!/usr/bin/env python

import argparse as ap, os, shutil as sh; from manim import *

Q2P = {
    "example" : "480p30",
    "fourk" : "2160p60",
    "high" : "1080p60",
    "low" : "480p15",
    "medium" : "720p30",
    "production" : "1440p60",
}

class Animation(Scene):
    def construct(self):
        axes = Axes(
            x_range=[-8, 8], y_range=[-1, 1, 10], axis_config={"include_tip": False}
        )
        labels = axes.get_axis_labels(x_label="t", y_label="\\psi(x)")

        sin_plot = axes.plot(lambda x: np.sin(x), color=BLUE)
        cos_plot = axes.plot(lambda x: np.cos(x), color=RED)

        self.play(Create(axes))
        self.play(Create(labels))
        self.play(Create(sin_plot))
        self.play(Transform(sin_plot, cos_plot))
        self.play(FadeOut(sin_plot))
        self.wait()

if __name__ == "__main__":
    # create the parser
    parser = ap.ArgumentParser(prog="Acorn Movie Renderer", description="Rendering script for the Quantum Acorn package.", add_help=False)

    # add the optional arguments
    parser.add_argument("-q", "--quality", type=str, default="production", help="The quality of the rendered image.")

    # parse arguments
    args = parser.parse_args()

    # throw error for invalid arguments
    if args.quality not in Q2P.keys():
        raise ValueError("INVALID QUALITY ARGUMENT")

    # set the manim options
    options = {
        "quality": args.quality + "_quality",
        "output_file" : "movie.mp4",
        "preview": False,
        "verbosity" : "ERROR"
    }

    # render the animation
    with tempconfig(options):
        Animation().render()

    # move the animation and remove the media directory
    sh.move("media/videos/{}/movie.mp4".format(Q2P[args.quality]), "movie.mp4"); sh.rmtree("media")
