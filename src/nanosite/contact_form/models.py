from django.db import models



class Contact_form(models.Model):
    user_name     = models.CharField(max_length=50, unique=True)
    email_address = models.EmailField(max_length=70)
    subject       = models.CharField(max_length=100)
    subject_type  = models.CharField(max_length=10, choices=(("SUGGESTION", "Suggestion"), ("FEEDBACK", "Feedback"), ("OTHER", "Other")),
                                     default="SUGGESTION")
    message       = models.TextField()