var photometry = [
{% for f in photometry %}
    {
        x: {{ photometry[f][0] }},
        y: {{ photometry[f][1] }},
        error_y: {
            type: 'data',
            array: {{ photometry[f][2] }},
            visible: true
        },
        name: '{{ f }}',
        text: {{ photometry[f][3] }},
        mode: 'markers',
        type: 'scatter',
{% if f[:2] == 's_' %}
        marker: {
            size: 10,
            symbol: 'circle',
        },
        opacity: .8,
{% else %}
        marker: {
            size: 10,
            symbol: 'square',
        },
        opacity: .5,
        visible: 'legendonly',
{% endif %}
    },
{% endfor %}
]
