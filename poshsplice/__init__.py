# -*- coding: utf-8 -*-

__author__ = 'Olga Botvinnik'
__email__ = 'olga.botvinnik@gmail.com'
__version__ = '0.1.0'


def wannabe():
    try:
        from IPython.display import YouTubeVideo
        YouTubeVideo("8X-2czaa3WA", autoplay=1, start=4)
    except:
        pass