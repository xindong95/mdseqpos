
class Filter:
    def __init__(self, fileName):
        self.f = open(fileName)
    
    def _nextLine(self):
        
        result = None
        
        while not result:
            result = self.f.next().strip()
            commentPos = result.find('#')
            if commentPos != - 1:
                result = result[:commentPos]

        return result
    
    def __iter__(self):
        return self
    
    def next(self):
        return self._nextLine()
