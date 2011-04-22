
# This is easy to solve with a simple tiny wrapper:
class Callable:
    def __init__(self, anycallable):
        self.__call__ = anycallable
                
class MyMath:
    def sum(a,b):
        if type(a)!=type(1):
            raise Exception
        if type(b)!=type(1):
            raise Exception
        return a+b
    def mul(a,b):
        pass
    def sub(a,b):
        pass
    #Making static methods
    sum = Callable(sum)
    mul = Callable(mul)
    sub = Callable(sub)
    
