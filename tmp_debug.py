import runpy, traceback
try:
    runpy.run_path('app/ui_gradio.py', run_name='__main__')
except Exception as exc:
    traceback.print_exc()
