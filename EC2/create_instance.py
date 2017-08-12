import boto3
import botocore
import time

ec2 = boto3.resource('ec2', region_name='us-east-1')
client = boto3.client('ec2')

# Create a security group
try:
    sg = ec2.create_security_group(GroupName='jupyter', Description='EC2 for Jupyter Notebook')
    response = client.authorize_security_group_ingress(GroupName='jupyter', IpPermissions=[{'PrefixListIds': [], 'UserIdGroupPairs': [], 'IpRanges': [{'CidrIp': '0.0.0.0/0'}], 'IpProtocol': 'tcp', 'Ipv6Ranges': [{'CidrIpv6': '::/0'}], 'ToPort': 8888, 'FromPort': 8888}])
    print("create a security group")
except botocore.exceptions.ClientError as e:
    sg = client.describe_security_groups(GroupNames=['jupyter'])
    print("the security group exist")

o = ec2.create_instances(ImageId='ami-622a0119', MinCount=1, MaxCount=1, InstanceType='i3.8xlarge', SecurityGroups=['jupyter'])
print_res = False
while (not print_res):
    time.sleep(1)
    for i in ec2.instances.filter(InstanceIds=[o[0].id]):
        if i.public_ip_address is not None:
            print("The public IP address: " + str(i.public_ip_address))
            print_res = True
