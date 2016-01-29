# -*- coding: utf-8 -*-

__author__ = 'Olga Botvinnik'
__email__ = 'olga.botvinnik@gmail.com'
__version__ = '0.1.0'


def wannabe():
    try:
        from IPython.display import YouTubeVideo, display
        vid = YouTubeVideo("8X-2czaa3WA", autoplay=1, start=4)
        display(vid)
    except ImportError:
        pass

__all__ = ['psi', 'events', 'gtf', 'junctions', 'psi', 'region', 'star',
           'util']
