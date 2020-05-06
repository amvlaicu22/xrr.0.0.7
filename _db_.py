import os
import sys
import inspect



############################################
from contextlib import contextmanager
from timeit import default_timer
@contextmanager
def Etimer():
    '''
    with Etimer() as elapsed:
        time.sleep(1)
        print(elapsed())
        time.sleep(1)
    print(elapsed())
    '''
    start = default_timer()
    elapser = lambda: default_timer() - start
    yield lambda: elapser()
    end = default_timer()
    elapser = lambda: end-start


def line(c="_", n=70): 
    return "".join([c for i in range(n)])
    
def linep(c="_", n=70): 
    print(line(c, n))
    
# ~ FRAM = lambda n=0: sys._getframe( 1+n )
# ~ LINE = lambda n=0: sys._getframe( 1+n ).f_lineno
# ~ CODE = lambda n=0: sys._getframe( 1+n ).f_code
# ~ NAME = lambda n=0: sys._getframe( 1+n ).f_code.co_name
# ~ FILE = lambda n=0: sys._getframe( 1+n ).f_code.co_filename
# ~ BASE = lambda f:   os.path.basename( f )
# ~ print(FRAM())
# ~ print(LINE())
# ~ print(CODE())
# ~ print(NAME())
# ~ print(FILE())
# ~ print(BASE(FILE()))

def insp_stack(cout, stk, *argv): 
    # *argv = [lev,[code, func, line, module]], ...
    # code = 4, func = 3, line = 2, module = 1 
    xx = []
    for arg in argv:
        lev = arg[0]
        wht = arg[1]
        for w in wht:
            x = stk[lev][w]
            if w == 1:  # module file
                x = os.path.basename(x)
            if w == 4:  # code
                x = x[0].lstrip().rstrip()
            # ~ if w == 3:  # func
                # ~ x = x
            xx.append(f"{x}")
    if cout:
        print( " | ".join(xx) )
    return xx

def insp(cout, *argv):  # d = dict(module=1, line=2, func=3, code=4)
    stk = inspect.stack()
    return insp_stack(cout, stk, *argv)
    
def insp_all(cout):
    stk = inspect.stack()
    x = []
    for lv in range(len(stk)):
        x += [f"{lv}"]
        x += insp_stack(cout, stk, [lv, [4,3,2,1]])
        x += ["\n"]
    return x


    

        
def _db_(*args, **kwargs):
    linep(c="~")
    print("!fpos:", end=''); insp(1, [2, [4, 3, 2, 1]])
    print("!call:", end=''); insp(1, [-1, [4, 3, 2, 1]])
    
    stk = inspect.stack()
    code = stk[1][4][0].lstrip().rstrip()   # _db_( v1, v2, ...)
    alist = code.lstrip('_db_(').rstrip(')').split(',')
    
    vlist = [*args] 
    if kwargs: vlist = vlist + [dict(kwargs)]
    if len(vlist):
        [ print(f"{a:>20} =  {v}") for a,v in zip(alist, vlist)]
    linep(c="~")
    


def argval(*args, **kwargs):
    args_repr = [repr(a) for a in args]                      # 1
    kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]  # 2
    sig = ", ".join(args_repr + kwargs_repr)           # 3
    return sig
    
def argval1(*args, **kwargs):
    s = ""
    for arg in args:  s += f"{arg}, "
    for key, value in kwargs.items():  s += f"({key}={value!r}) , "
    return s


def argval2(*args, **kwargs):
    return f"{[*args] + [dict(kwargs)]}"

def funcsig(func, *args, **kwargs):
    return f"{func.__name__!r}({argval(*args, **kwargs)})"

# ~ decorator template
import functools
def Decorator(func, func1=None, func2=None):
    # ~ fdec = insp(1, [1,[3]])[0] 
    # ~ call = insp(1, [-1, [4, 3, 2, 1]]) 
    fpos = insp(0, [2, [4, 3, 2, 1]])
    pos  = " | ".join(fpos)
    
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        fsig = f"{func.__name__!r}({argval(*args, **kwargs)})"
        # ~ if func1: fsig += f" | {func1.__name__!r}"
        # ~ if func2: fsig += f" | {func2.__name__!r}"
        sig = " | ".join( [f"{pos}", fsig  ] )
        print(f"{sig}")
        if func1: func1()
        value = func(*args, **kwargs)
        if func2: func2()
        print(f"{sig}")
        linep()
        return value
    return wrapper_decorator

import functools
import time
def Dtimer(func, *args0, **kwargs0):
    # ~ fdec = insp(0, [1,[3]])[0] 
    fpos0 = "|".join( insp(0, [ 2, [4, 2, 1]]) )
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        linep("-")
        fpos1 = "|".join( insp(0, [ 2, [4, 2, 1]]) )
        # ~ fsig = f"{func.__name__!r}({argval(*args, **kwargs)})"
        sig = f"{fpos1} :: {fpos0}"
        print(f"{sig}")    
        start_time = time.perf_counter()    # 1
        #############################
        value = func(*args, **kwargs)
        #############################
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"dT ::= {run_time:.4f} s")
        print(f"{sig}")
        linep("-")
        return value
    return wrapper_timer






class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    """
    usage: 
    1) with Timer("T1: with Timer() ...") as test1: { code block}
    2) @Timer()
       def func(): ...
    3) t = Timer(text="ABC"); t.start(); ...
    * possible changes to __init__
    __init__(self, text="", logger=print, format="Timer {} elapsed time: {:0.4f} seconds")
        self.log = logger
        self.frm = format
        if self.log:
            self.log(self.form.format(self.name, elapsed_time))
    """
    
    timers = dict()
    name = None

    def _insp(self):  
        insp(0, [2,[3]], [-1,[4,3,2,1]] )
    
    def __init__(self, text=""):
        linep()
        self._insp()
        self._start_time = None
        self.text = text
        name = [f">Timer"]
        name += insp(0, [2, [4, 3, 2, 1]]) 
        if not self.name :
            self.name = " | ".join(name)
        self.timers.setdefault(self.name, 0)
        # ~ print(self.name)
        # ~ Add new named timers to dictionary of timers
        # ~ only sets the value if name is not already defined in the dictionary. 
        # ~ If name is already used in .timers , then the value is left untouched
        # ~ insp_all(1)
            
    def start(self):
        self._insp()
        print(self.name)
        
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")
        self._start_time = time.perf_counter()
        
    def stop(self):
        self._insp()
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")
        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        
        if self.name:
            self.timers[self.name] += elapsed_time
        
        print(f"elapsed_time::= {elapsed_time:.4f} s")
        print(f"{self.name}")
        linep()
        return elapsed_time
        
    
        
    def __call__(self, func):
        self._insp()
        """
        make the class callable in order to be used as decorator
        called when using Timer as a decorator
        """
        @functools.wraps(func)
        def wrapper_timer(*args, **kwargs):
            with self:
                return func(*args, **kwargs)
        return wrapper_timer
        
    def __enter__(self):
        self._insp()
        """ called when using: with Timer() as t: """
        self.start()
        return self
    
    def __exit__(self, *exc_info):
        self._insp()
        '''
        def __exit__(self, exc_type, exc_value, exc_tb):
        =================================================
        __exit__() takes three arguments: exc_type , exc_value , and exc_tb .
        These are used for error handling within the context manager, and mirror the
        return values of sys.exc_info() . If an exception happens while the block is being
        executed, then your code calls .__exit__() with the type of the exception, an
        exception instance, and a traceback object. Often, you can ignore these in your
        context manager, in which case .__exit__() is called before the exception is
        reraised:
        '''
        self.stop()
    
    
    
from contextlib import ContextDecorator
class TimerCD(Timer, ContextDecorator):
    '''
    This change is just syntactic sugar for any construct of the following form:
    def f():
        with TimerCD():
            # Do stuff
            
    ContextDecorator lets you instead write:
    
    @TimerCD()
    def f():
        # Do stuff
    '''
    def __init__(self, *args, **kwargs):
        
        self.name = None
        name = [f">TimerCDs"]
        # ~ name += insp(0, [-1, [4, 3, 2, 1]]) 
        name += insp(0, [2, [4, 3, 2, 1]]) 
        self.name = " | ".join(name)
        self.timers.setdefault(self.name, 0)
        # ~ print(self.name)
        super().__init__(*args, **kwargs)
    



'''
# ~ from icecream import ic
# ~ import q
# ~ import traceback
# ~ from numba import jit
# ~ import snoop
# ~ from icecream import ic
# ~ from birdseye import eye
# ~ import peepshow
# ~ from decorators import debug
# ~ @debug
#############################
# ~ import pdb; 
# ~ pdb.set_trace()
# ~ answ = input('Enter pdb? (Ctrl-c / <Enter> / y=pdb):')         #py3
# ~ print(f"your choice is {answ.lower()}")
# ~ if "y" == answ.lower() : print('entering pdb.set_trace()')
###################################
# ~ import inspect
# ~ def foo(arg1,arg2):         
    # ~ #do something with args 
    # ~ a = arg1 + arg2         
    # ~ return a  
# ~ lines = inspect.getsource(foo)
# ~ print(lines)
######################################
# ~ import dis
# ~ dis.dis(foo)
######################################
# ~ from dill.source import getsource

# ~ def add(x,y):
   # ~ return x+y
# ~ print getsource(add)

# ~ squared = lambda x:x**2
# ~ print getsource(squared)

# ~ class Foo(object):
   # ~ def bar(self, x):
     # ~ return x*x+x
 
# ~ f = Foo()
# ~ print getsource(f.bar)
##############################################
# ~ f = lambda x: x.method_foo
# ~ in  py3 
# ~ The func_code attribute of functions has been renamed to __code__, 
# ~ ff = f.__name__   
# ~ ff = f.__code__.co_names        # ~ print(ff)   # ~ ('method_foo',)
# ~ ff = f.__code__.co_names[0]     # ~ print(ff)   # ~ 'method_foo'
# ~ ff = f.__code__.co_filename
# ~ ff = f.__code__.co_firstlineno

# ~ in  py2
# ~ ff=f.func_name  
# ~ ff=f.func_code.co_names in py2
# ~ ff=f.func_code.co_names[0] in py2
# ~ print ff
################################################
# ~ f formatting
# ~ ----------------
# ~ print(d, end='', flush=True)
# ~ f ' <text> { <expression> <optional !s, !r, or !a> <optional : format specifier> } <text> ... '
# ~ Fahrenheit = [32, 60, 102]
# ~ [f'{((x - 32) * (5/9)):.2f} Celsius' for x in Fahrenheit]
# ~ >> ['0.00 Celsius', '15.56 Celsius', '38.89 Celsius']



'''


