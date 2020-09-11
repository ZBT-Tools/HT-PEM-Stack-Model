from . import widget_set
from . import button


class WidgetFactory:

    def __init__(self):
        self._widget_set_factory = widget_set.WidgetSetFactory()
        self._button_factory = button.ButtonFactory()

    def create(self, frame, **kwargs):
        widget_type = kwargs.get('type', None)
        if 'Set' in widget_type:
            return self._widget_set_factory.create(frame, **kwargs)
        elif 'Button' in widget_type and not 'Set' in widget_type:
            return self._button_factory.create(frame, **kwargs)
        elif widget_type == 'Label':
            kwargs.pop('type', None)
            return widget_set.Label(frame, **kwargs)
        else:
            raise NotImplementedError('type of widget not implemented')
