from django.shortcuts import render
from django.http import HttpResponse


# Website views
def home_view(request, *args, **kwargs):
    return render(request, "home.html", {})

def contact_view(request, *args, **kwargs):
    return render(request, "contact.html", {})

def siteStats_view(request, *args, **kwargs):
    return render(request, "stats.html", {})

def howToCite_view(request, *args, **kwargs):
    return render(request, "cite.html", {})

def report_view(request, *args, **kwargs):
    return render(request, "report.html", {})

def overview_view(request, *args, **kwargs):
    return render(request, "overview.html", {})

def howToUse_view(request, *args, **kwargs):
    return render(request, "usage.html", {})