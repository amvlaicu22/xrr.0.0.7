import _db_
from _db_ import _db_
from _db_ import *


def elapsed_timer_test():
    with Etimer() as elapsed:
        print("test elapsed_timer")
        print("\twith elapsed_timer() as elapsed:\n\t ... print(elapsed())")
        time.sleep(0.1)
        print(elapsed())
        time.sleep(0.1)
        print(elapsed())
        time.sleep(0.1)
    print(elapsed())
    linep()
    
def _db_test():
    _db_()
    def _db_test2():
        a = 1
        _db_("inside")
        _db_('test', a, 1, 2, 3)
    _db_test2()
    
def deco_test():
    @Decorator
    def decorator_test():
        print("inside decorator_test")
        
    decorator_test()
    
    #####################
    __t0 = time.time()
    def ff(a=1): 
        f = insp(1, [1,[3]])      
        # ~ print(istack_func())
        time.sleep(0.1)
        
    def f1():   
        f = insp(1, [1,[3]])
        # ~ global __t0
        __t0 = time.time()
        
    def f2():   
        f = insp(1, [1,[3]])
        # ~ global __t0
        print(time.time()-__t0)
        
    dt = Decorator(ff, f1, f2)
    dt(a=2)
    
    
def timer_deco_test():
    
    @Dtimer
    def timer_test():
        print("inside timer_test")
        
    timer_test()
    
    #######################
    import math
    
    
    @Dtimer
    def approximate_e(terms=100):
        fact = math.factorial
        return sum(1 / fact(n) for n in range(terms))
    
    approximate_e(5)
    approximate_e(500)
    
    #################################
    
    
    print("# Apply a decorator to a standard library function")
    print("# fact = Dtimer(math.factorial)")
    print("# print(fact(3))")
    # ~ linep()
    fact = Dtimer(math.factorial)
    print(fact(3))
    
    print("# return to standard")
    print("# fact = math.factorial")
    print("# print(fact(3))")
    # ~ linep()
    fact = math.factorial
    print(fact(3))
    # ~ linep()
    
    
def Timer_class_test():
    @Timer("deco")
    def waste_some_time(num_times):
        for _ in range(num_times):
            sum([i**2 for i in range(10000)])
            
    waste_some_time(10)

    #################################
    
    
    def with_Timer_test():
        with Timer("with") as t:
            for i in range(10):
                time.sleep(0.01)

    with_Timer_test()
    
    ############################
    t = Timer("inst")
    t.start()
    time.sleep(0.1)
    t.stop()
    
    ######################

    with Timer("with") as t:
        time.sleep(0.1)

    @TimerCD(" _db_ ")
    def testfunc_Timer_deco_with_db():
        m = "abc"
        _db_(1,m,len("abs"), 3)
        
    testfunc_Timer_deco_with_db()
    #############################
    
    
if __name__ == "__main__":
        
    # ~ elapsed_timer_test()
    
    # ~ _db_test()

    # ~ deco_test()
    
    timer_deco_test()
    
    # ~ Timer_class_test()
    
    
    
    
