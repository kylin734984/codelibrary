import argparse
'''
module API demo
'''
import logging
import sys
import os
from functools import wraps
import time

class options():
    def __int__(self):
        self.infp = None
        self.outfp = None

class Logger():
    def __init__(self, name = None):
        self.logger = logging.getLogger(name) if name else logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.formatter = logging.Formatter('%(message)s')
        self.f_handler = logging.FileHandler('{0}.log'.format(__name__))
        self.f_handler.setFormatter(self.formatter)
        self.f_handler.setLevel(logging.INFO)
        self.logger.addHandler(self.f_handler)
        self.console_handler = logging.StreamHandler()
        self.console_handler.setFormatter(self.formatter)
        self.console_handler.setLevel(logging.INFO)
        self.logger.addHandler(self.console_handler)
        self.logger.setLevel(logging.INFO)

    def rename(self, outprefix):
        self.f_handler.flush()
        self.f_handler.close()
        self.logger.removeHandler(self.f_handler)
        if isinstance(outprefix,str):
            try:
                os.remove('{0}.log'.format(outprefix))
            except:
                pass
            try:
                os.rename('{0}.log'.format(__name__), '{0}.log'.format(outprefix))
            except:
                print('Can\'t rename log file: {0}'.format('{0}.log'.format(__name__)))



def timeit(logger):
    def decorator(fun):
        @wraps(fun)  #为了保留被装饰函数的函数名和帮助文档信息
        def wrapper(*args,**kwargs):
            """这是一个wrapper函数"""
            start_time = time.time()
            logger.logger.info('Start time: {0}'.format(time.strftime('%H:%M:%S', time.localtime(start_time))))
            res = fun(*args,**kwargs)
            end_time = time.time()
            use = time.gmtime(end_time - start_time)
            logger.logger.info('Complete time: {0}\nEclipsed time: {1}'.format(time.strftime('%H:%M:%S', time.localtime(end_time)),
                                                                  time.strftime('%H:%M:%S', use)))
            return res
        return wrapper
    return decorator

def parse_arg():
    parser = argparse.ArgumentParser(description='A module API demo')
    parser.add_argument('-i', '--input', dest='infp', required=True, default=None, metavar='FILE', help='Input file')
    parser.add_argument('-o', '--output', dest='outfp', default=None, metavar='FILE', help='Output file')
    return parser.parse_args()

def methodbody(opts, logger):
    ''':param
    Mian method
    '''
    logger.logger.info('This is a piece of information about logging module')


logger = Logger()

@timeit(logger)
def init(opts, logger):
    methodbody(opts, logger)
    logger.rename(opts.outfp)

if __name__ == '__main__':
    if len(sys.argv)==1:
        opts = options()
        opts.infp = 'a'
        opts.outfp = 'b'
        init(opts, logger)
    else:
        opts = parse_arg()
        init(opts, logger)



