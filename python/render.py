#!/usr/bin/env python

import argparse, os, shutil; from manim import *

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
        circle = Circle()
        circle.set_fill(PINK, opacity=0.5)
        self.play(Create(circle))

if __name__ == "__main__":
    # create the parser
    parser = argparse.ArgumentParser(prog="Acorn Movie Renderer", description="Rendering script for the Quantum Acorn package.", add_help=False)

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
        "preview": False,
        "verbosity" : "ERROR"
    }

    # render the animation
    with tempconfig(options):
        Animation().render()

    # move the animation and remove the media directory
    shutil.move("media/videos/{}/Animation.mp4".format(Q2P[args.quality]), "movie.mp4"); shutil.rmtree("media")
