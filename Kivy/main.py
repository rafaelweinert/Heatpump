import kivy
from kivy.app import App
from kivy.uix.label import Label

from kivy.lang import Builder

# this restrict the kivy version i.e

from kivy.uix.boxlayout import BoxLayout


###############################################

# creating the root widget used in .kv file
class KVBL(BoxLayout):
    '''
        no need to do anything here as
        we are building things in .kv file
    '''
    pass


#################################################
# class in which name .kv file must be named KVBoxLayout.kv.
class KVBoxLayoutApp(App):

    def build(self):
        # returning the instance of KVBL class
        return KVBL()


##################################################

# creating the object root for BoxLayoutApp() class
root = KVBoxLayoutApp()

# run function runs the whole program
# i.e run() method which calls the
# target function passed to the constructor.
root.run()

